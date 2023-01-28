from typing import Dict, List, Union

import datetime
import inspect
import tempfile
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import dicom2nifti
import highdicom as hd
import nibabel as nb
import numpy as np
import pandas as pd
import pkg_resources
import pydicom
import yaml
from colormath.color_conversions import convert_color
from colormath.color_objects import LabColor, sRGBColor
from fhir.resources import construct_fhir_element
from highdicom.sr import content, templates
from nilearn import image
from nilearn.image import resample_to_img
from pydicom.sr.codedict import codes
from pydicom.uid import generate_uid

TEMPLATE_PATH = pkg_resources.resource_filename("pinr", "templates")
template_file = Path(TEMPLATE_PATH) / "fs-aseg.yml"

CMET_ROOT_UID = "1.2.826.0.1.3680043.10.538."

# #TODO: update names/versions?
SEG_DEVICE_NAME = "FreeSurfer"
DEVICE_MANUFACTURER = "CorticoMetrics"
DEVICE_OBSERVER_UID = generate_uid(CMET_ROOT_UID)
SR_DEVICE_NAME = "FreeSurfer processed by pinr"
SR_MANUFACTURER = "CorticoMetrics"
SR_VERSION = "0.0.1"  # #TODO import __version__

FS_ALGORITHM = "FreeSurfer"
FS_VERSION = "7"
SERIAL_NUMBER = f"pinr-{SR_VERSION}"

REPORT_NAME = "FreeSurfer ASeg Measurement Report, created by pinr"

observation_datetime = datetime.datetime.now()
patient_gender_map = dict(M="male", F="female", O="other")

# #TODO: clean up and make less FS specific
# FreeSurfer has a weird format for cortex and white mater, need to map these to be standard
wholebrain_structure_id = {
    "lhCerebralWhiteMatterVol": 2,
    "lhCortexVol": 3,
    "rhCerebralWhiteMatterVol": 41,
    "rhCortexVol": 42,
}

wholebrain_structure_name = {
    "lhCerebralWhiteMatterVol": "Left-Cerebral-White-Matter",
    "lhCortexVol": "Left-Cerebral-Cortex",
    "rhCerebralWhiteMatterVol": "Right-Cerebral-White-Matter",
    "rhCortexVol": "Right-Cerebral-Cortex",
}


@dataclass
class FreeSurferSeg:
    """Stores Volumetric Data from an aseg.stats file.

    Attributes
    ----------
    t1w_dicom_file : Union[str, Path]
        One T1w DICOM file, part of the series used as input to FreeSurfer
    aseg_file : Union[str, Path]
        aseg.mgz Segmentation file created by FreeSurfer
    aseg_stats_file : Union[str, Path]
        aseg.stats file output by FreeSurfer.
    dicom_seg_output_file : Union[str, Path, None]
        output file name for DICOM Seg
    dicom_sr_output_file : Union[str, Path, None]
        output file name for DICOM SR
    fhir_output_file : Union[str, Path, None]
        output file name for FHIR Diagnostic Report JSON
    dicom_seg_schema_file : Union[str, Path, None]
        Schema file used to provided coded inputs for each structure
    """

    t1w_dicom_file: Union[str, Path]
    aseg_file: Union[str, Path]
    aseg_stats_file: Union[str, Path]

    dicom_seg_output_file: Union[str, Path, None] = None
    dicom_sr_output_file: Union[str, Path, None] = None
    fhir_output_file: Union[str, Path, None] = None

    dicom_seg_schema_file: Union[str, Path, None] = None

    @staticmethod
    def _get_seg_attr_vals(seg_attr, ordered=True):
        if ordered:
            num = seg_attr["ordered_label_id"]
        else:
            num = seg_attr["label_id"]
        label = seg_attr["segment_description"].replace("-", " ").replace("_", " ")
        property_type_meaning = (
            seg_attr["segmented_property_type"]
            .title()
            .replace(" ", "")
            .replace("_", "")
        )
        try:
            # Sometimes there are multiple codes for the same structure
            # Try using the "CID 7140: Brain Structure for Volumetric Measurement" code first
            # https://dicom.nema.org/medical/dicom/current/output/chtml/part16/sect_cid_7140.html
            property_type_code = getattr(codes.cid7140, property_type_meaning)
        except AttributeError:
            property_type_code = getattr(getattr(codes, "SCT"), property_type_meaning)
        # #NOTE: unsure on the best laterality
        # FHIR BodyStructure may require "Unilateral Left/Right", but there's no code for UnilateralRight in pydicom
        # https://www.hl7.org/fhir/valueset-bodystructure-relative-location.html
        if label.startswith("Left"):
            laterality = codes.SCT.Left
        elif label.startswith("Right"):
            laterality = codes.SCT.Right
        else:
            laterality = None
        structure_uid = seg_attr["structure_tracking_uid"]
        measure_uid = seg_attr["measure_tracking_uid"]
        return num, label, property_type_code, laterality, structure_uid, measure_uid

    @cached_property
    def volume_measurements(self) -> Dict[str, float]:
        """Load data from a FreeSurfer aseg.stats file."""
        # #TODO: clean up and make less FS specific?
        column_headers = [
            "SegId",
            "NVoxels",
            "Volume_mm3",
            "StructName",
            "normMean",
            "normStdDev",
            "normMin",
            "normMax",
            "normRange",
        ]

        df = pd.read_table(
            self.aseg_stats_file,
            delim_whitespace=True,
            header=None,
            comment="#",
            index_col=0,
            names=column_headers,
        )
        # #HACK: read in wholebrain wm and cortex
        with open(self.aseg_stats_file, "r", encoding="utf-8") as f:
            cortical_data = []
            for x in f:
                if x.startswith("# Measure "):
                    structure_name = x.split(",")[1].strip()
                    if structure_name in wholebrain_structure_name:
                        cortical_data.append(
                            {
                                "SegId": wholebrain_structure_id[structure_name],
                                "NVoxels": np.nan,
                                "Volume_mm3": float(x.split(",")[-2]),
                                "StructName": wholebrain_structure_name[structure_name],
                                "normMean": np.nan,
                                "normStdDev": np.nan,
                                "normMin": np.nan,
                                "normMax": np.nan,
                                "normRange": np.nan,
                            }
                        )
        cortical_df = pd.DataFrame(cortical_data)
        df = pd.concat([df, cortical_df]).reset_index(drop=True)

        # Only include structures in our schema:
        data = []
        for x in self._seg_attrs:
            seg = x["segment_description"]
            try:
                values = df.loc[df.StructName == seg][
                    ["StructName", "Volume_mm3"]
                ].values[0]
            except IndexError:
                # #HACK: Older versions of FS have Thalmus-Proper as struct name
                if seg.endswith("Thalamus"):
                    vol = df.loc[df.StructName == f"{seg}-Proper"]["Volume_mm3"].values[
                        0
                    ]
                    values = [seg, vol]
            data.append(values)
        return dict(data)

    @cached_property
    def t1w_image_datasets(self) -> List[pydicom.FileDataset]:
        dcm_dir = Path(self.t1w_dicom_file).parent
        dcm = pydicom.dcmread(str(self.t1w_dicom_file))
        dcm_dict = {}
        for x in dcm_dir.glob("*"):
            try:
                d = pydicom.dcmread(x)
                if d.SeriesInstanceUID == dcm.SeriesInstanceUID:
                    try:
                        _ = d.SpacingBetweenSlices
                    except AttributeError:
                        d.SpacingBetweenSlices = d.SliceThickness
                    dcm_dict[int(d.InstanceNumber)] = d
            except pydicom.errors.InvalidDicomError:
                pass
        sorted_keys = sorted(dcm_dict)
        dcm_data = []
        for y in sorted_keys:
            dcm_data.append(dcm_dict[y])
        return dcm_data

    @cached_property
    def t1w_image_dataset(self) -> pydicom.FileDataset:
        return pydicom.dcmread(str(self.t1w_dicom_file))

    @cached_property
    def resampled_aseg(self) -> nb.Nifti1Image:
        aseg_image = image.load_img(str(self.aseg_file))
        aseg_image = nb.Nifti1Image(
            aseg_image.dataobj, aseg_image.affine, aseg_image.header
        )
        with tempfile.TemporaryDirectory() as _temp_dir:
            t1w_nii = dicom2nifti.convert_dicom.dicom_array_to_nifti(
                self.t1w_image_datasets,
                Path(_temp_dir, "t1w.nii.gz"),
                reorient_nifti=False,
            )
            t1w_image = image.load_img(str(t1w_nii["NII_FILE"]))
            resampled = resample_to_img(aseg_image, t1w_image, interpolation="nearest")

            # use ordered_label_ids to change FS label numbers to sequential
            for x in self._ordered_seg_attrs:
                label = x["label_id"]
                ordered_label = x["ordered_label_id"]
                resampled_data = resampled.dataobj
                resampled_data[resampled_data == label] = ordered_label
                resampled = nb.Nifti1Image(
                    resampled_data, resampled.affine, resampled.header
                )
            return resampled

    @cached_property
    def segmentation_schema(self) -> dict:
        if self.dicom_seg_schema_file is None:
            schema_file = template_file
        else:
            schema_file = self.dicom_seg_schema_file
        with open(schema_file, "r", encoding="utf-8") as f:
            aseg_schema = yaml.safe_load(f)

        seg_attrs = aseg_schema["segment_attributes"].copy()
        for i, x in enumerate(seg_attrs):
            # add in recommended color, in CIELab format
            rgb = sRGBColor(*x["recommended_rgb_value"], is_upscaled=True)
            lab = convert_color(rgb, LabColor)
            # pylint: disable=no-value-for-parameter
            CIELabColor = hd.color.CIELabColor(*lab.get_value_tuple())
            aseg_schema["segment_attributes"][i][
                "recommended_cielab_color"
            ] = CIELabColor
            # Use the same UID for the Seg and SR structure
            aseg_schema["segment_attributes"][i][
                "structure_tracking_uid"
            ] = generate_uid(CMET_ROOT_UID)
            # Use the same UID for the measurement
            aseg_schema["segment_attributes"][i]["measure_tracking_uid"] = generate_uid(
                CMET_ROOT_UID
            )
        return aseg_schema

    @cached_property
    def seg(self):
        descriptions = []
        for x in self._ordered_seg_attrs:
            num, label, code, _, structure_uid, _ = self._get_seg_attr_vals(x)

            d = hd.seg.SegmentDescription(
                segment_number=num,
                segment_label=label,
                segmented_property_category=hd.sr.CodedConcept(
                    "91723000", "SCT", "Anatomical Structure"
                ),
                segmented_property_type=code,
                algorithm_type=hd.seg.SegmentAlgorithmTypeValues.AUTOMATIC,
                algorithm_identification=hd.AlgorithmIdentificationSequence(
                    name=FS_ALGORITHM,
                    version=FS_VERSION,
                    family=codes.cid7162.ArtificialIntelligence,
                ),
                tracking_uid=structure_uid,
                tracking_id=label,
            )
            descriptions.append(d)
        seg = hd.seg.Segmentation(
            source_images=self.t1w_image_datasets,
            pixel_array=np.swapaxes(self.resampled_aseg.dataobj.astype("uint16"), 2, 0),
            segmentation_type=hd.seg.SegmentationTypeValues.BINARY,
            segment_descriptions=descriptions,
            series_instance_uid=generate_uid(CMET_ROOT_UID),
            series_number=777,
            sop_instance_uid=generate_uid(CMET_ROOT_UID),
            instance_number=1,
            manufacturer=DEVICE_MANUFACTURER,
            manufacturer_model_name=SEG_DEVICE_NAME,
            software_versions=FS_VERSION,
            device_serial_number=SERIAL_NUMBER,
            omit_empty_frames=True,
            content_description="FreeSurfer ASeg DICOM Segmentation, created pinr",
            series_description="Segmentation",
        )
        seg_sequence = seg.SegmentSequence
        # add in the recommended color values
        for i, x in enumerate(self._ordered_seg_attrs):
            seg_sequence[i].RecommendedDisplayCIELabValue = list(
                x["recommended_cielab_color"].value
            )
        if self.dicom_seg_output_file is not None:
            seg.save_as(self.dicom_seg_output_file)
        return seg

    @cached_property
    def _seg_attrs(self) -> List[dict]:
        return self.segmentation_schema["segment_attributes"].copy()

    @cached_property
    def _ordered_seg_attrs(self):
        ordered_seg_attrs = []
        ordered_labels = self._ordered_label_map
        for x in self._seg_attrs:
            label = x["label_id"]
            if label in ordered_labels:
                x["ordered_label_id"] = ordered_labels[label]
                ordered_seg_attrs.append(x)
        return ordered_seg_attrs

    @cached_property
    def _ordered_label_map(self):
        # get ordered label IDs, not just FS label IDs, only for vols > 0
        vols = self.volume_measurements
        label_ids = sorted(
            [
                x["label_id"]
                for x in self._seg_attrs
                if vols[x["segment_description"]] > 0
            ]
        )
        ordered_ids = range(1, len(label_ids) + 1)
        label_map = dict(zip(label_ids, ordered_ids))
        return label_map

    @cached_property
    def sr(self):
        observer_device_context = hd.sr.ObserverContext(
            observer_type=codes.DCM.Device,
            observer_identifying_attributes=hd.sr.DeviceObserverIdentifyingAttributes(
                uid=DEVICE_OBSERVER_UID,
                name=SR_DEVICE_NAME,
                manufacturer_name=DEVICE_MANUFACTURER,
            ),
        )
        observation_context = hd.sr.ObservationContext(
            observer_device_context=observer_device_context,
        )
        imaging_measurements = []
        vols = self.volume_measurements
        for x in self._ordered_seg_attrs:
            (
                num,
                label,
                code,
                laterality,
                structure_uid,
                measure_uid,
            ) = self._get_seg_attr_vals(x)
            tid300_measurement = [
                templates.Measurement(
                    name=codes.SCT.Volume,
                    value=vols[x["segment_description"]],
                    unit=codes.UCUM.CubicMillimeter,
                    tracking_identifier=templates.TrackingIdentifier(
                        uid=measure_uid,
                        identifier=f"{label} Volume",
                    ),
                    # #TODO: add properties for normative measures
                    # referenced_images = [
                    # hd.sr.SourceImageForMeasurement.from_source_image(y) for y in self.t1w_image_datasets
                    # ]
                )
            ]
            tid1411_vol_roi = hd.sr.VolumetricROIMeasurementsAndQualitativeEvaluations(
                tracking_identifier=templates.TrackingIdentifier(
                    uid=structure_uid,
                    identifier=label,
                ),
                referenced_segment=hd.sr.ReferencedSegment.from_segmentation(
                    self.seg, num
                ),
                finding_type=hd.sr.CodedConcept(
                    "91723000", "SCT", "Anatomical Structure"
                ),
                algorithm_id=hd.sr.AlgorithmIdentification(
                    name=FS_ALGORITHM, version=FS_VERSION
                ),
                finding_sites=[
                    content.FindingSite(
                        anatomic_location=code,
                        laterality=laterality,
                    )
                ],
                measurements=tid300_measurement,
                # #TODO: add time_point_context for longitidunal
            )

            imaging_measurements.append(tid1411_vol_roi)
        measurement_report = templates.MeasurementReport(
            observation_context=observation_context,
            procedure_reported=codes.SCT.ImagingProcedure,
            imaging_measurements=imaging_measurements,
            title=codes.DCM.ImagingMeasurementReport,
            referenced_images=self.t1w_image_datasets,
        )
        sr_dataset = hd.sr.Comprehensive3DSR(
            evidence=self.t1w_image_datasets + [self.seg],
            content=measurement_report[0],
            series_number=888,
            series_instance_uid=generate_uid(CMET_ROOT_UID),
            sop_instance_uid=generate_uid(CMET_ROOT_UID),
            instance_number=1,
            manufacturer=SR_MANUFACTURER,
            is_complete=True,
            is_verified=False,
            institution_name=SR_MANUFACTURER,
            series_description="Volumetric Measurements",
        )

        if self.dicom_sr_output_file is not None:
            sr_dataset.save_as(self.dicom_sr_output_file)
        return sr_dataset

    @cached_property
    def _device(self):
        """See here for the DICOM SR specific requirements.

        https://build.fhir.org/ig/HL7/dicom-sr/StructureDefinition-tid-4019-algorithm-identification.html
        """
        json_data = {
            "id": "SRDevice",
            "resourceType": "Device",
            "identifier": [
                {
                    "type": {
                        "coding": [
                            {
                                "code": "121012",
                                "display": "Device Observer UID",
                                "system": "http://dicom.nema.org/resources/ontology/DCM",
                            }
                        ]
                    },
                    "use": "usual",
                    "value": DEVICE_OBSERVER_UID,
                },
                {
                    "type": {
                        "coding": [
                            {
                                "code": "121013",
                                "display": "Device Observer Name",
                                "system": "http://dicom.nema.org/resources/ontology/DCM",
                            }
                        ]
                    },
                    "use": "usual",
                    "value": SR_DEVICE_NAME,
                },
                {
                    "type": {
                        "coding": [
                            {
                                "code": "121014",
                                "display": "Device Observer Manufacturer",
                                "system": "http://dicom.nema.org/resources/ontology/DCM",
                            }
                        ]
                    },
                    "use": "usual",
                    "value": SR_MANUFACTURER,
                },
            ],
            "deviceName": [{"name": SR_DEVICE_NAME, "type": "model-name"}],
            "version": [{"value": SR_VERSION}],
        }
        element = construct_fhir_element("Device", json_data)
        return element

    @cached_property
    def _algorithm_identification(self):
        """See here for the DICOM SR specific requirements.

        https://build.fhir.org/ig/HL7/dicom-sr/StructureDefinition-tid-4019-algorithm-identification.html
        """
        json_data = {
            "id": "AlgorithmIdentification",
            "resourceType": "Device",
            "deviceName": [{"name": FS_ALGORITHM, "type": "model-name"}],
            "version": [{"value": FS_VERSION}],
        }
        element = construct_fhir_element("Device", json_data)
        return element

    @cached_property
    def _patient(self):
        dcm_details = self._t1w_image_properties
        json_data = {
            "id": "Patient",
            "resourceType": "Patient",
            "identifier": [
                {
                    "type": {
                        "coding": [
                            {
                                "system": "http://terminology.hl7.org/CodeSystem/v2-0203",
                                "code": "ACSN",
                            }
                        ]
                    },
                    "use": "usual",
                    "value": dcm_details["acsn"],
                }
            ],
            "name": [{"text": str(dcm_details["PatientName"])}],
            "birthDate": pydicom.valuerep.DA(dcm_details["PatientBirthDate"]),
            "gender": patient_gender_map[dcm_details["PatientSex"]],
        }
        element = construct_fhir_element("Patient", json_data)
        return element

    @cached_property
    def _imaging_study(self):
        patient_ref = f"#{self._patient.id}"
        dcm_details = self._t1w_image_properties
        json_data = {
            "id": "ImagingStudy-T1w",
            "resourceType": "ImagingStudy",
            "identifier": [
                {
                    "type": {
                        "coding": [
                            {
                                "code": "study-instance-uid",
                                "display": "Study Instance UID",
                                "system": "http://hl7.org/fhir/uv/dicom-sr/CodeSystem/dicom-identifier-type",
                            }
                        ]
                    },
                    "use": "usual",
                    "value": dcm_details["StudyInstanceUID"],
                },
                {
                    "type": {
                        "coding": [
                            {
                                "code": "110180",
                                "display": "Study Instance UID",
                                "system": "http://dicom.nema.org/resources/ontology/DCM",
                            }
                        ]
                    },
                    "use": "usual",
                    "value": dcm_details["StudyInstanceUID"],
                },
                {
                    "type": {
                        "coding": [
                            {
                                "code": "ACSN",
                                "system": "http://terminology.hl7.org/CodeSystem/v2-0203",
                            }
                        ]
                    },
                    "use": "usual",
                    "value": dcm_details["acsn"],
                },
            ],
            "status": "available",
            "subject": {"reference": patient_ref},
            "procedureCode": [
                {
                    "coding": [
                        {
                            "code": "30657-1",
                            "display": "Brain",
                            "system": "MR Brain WO contr",
                        },
                    ],
                    "text": "MR Brain WO contrast",
                }
            ],
            "started": dcm_details["study_datetime"],
            "series": [
                {
                    "uid": dcm_details["SeriesInstanceUID"],
                    "number": dcm_details["SeriesNumber"],
                    "numberOfInstances": dcm_details["num_instances"],
                    "modality": {
                        "code": "MR",
                        "display": "Magnetic Resonance",
                        "system": "http://dicom.nema.org/resources/ontology/DCM",
                    },
                    "bodySite": {
                        "code": "RID6434",
                        "display": "Brain",
                        "system": "http://www.radlex.org",
                    },
                },
            ],
        }
        element = construct_fhir_element("ImagingStudy", json_data)
        return element

    @cached_property
    def _body_structures(self):
        # NOTE: the spec uses an updated version of BodyStructure, with includedStructure instead of location
        # https://build.fhir.org/ig/HL7/dicom-sr/StructureDefinition-dicom-sr-finding-site-body-structure.html
        # http://build.fhir.org/bodystructure.html
        # #TODO:
        patient_ref = f"#{self._patient.id}"
        body_structures = []
        for x in self._ordered_seg_attrs:
            (_, label, code, laterality, structure_uid, _) = self._get_seg_attr_vals(x)
            if code.scheme_designator == "SCT":
                code_system = "http://snomed.info/sct"
            elif code.scheme_designator == "DCM":
                code_system = "http://dicom.nema.org/resources/ontology/DCM"
            else:
                code_system = None  # #TODO: handle different systems
            stucture_id = f"{label.title().replace(' ', '').replace('_', '')}-Structure"
            json_data = {
                "id": stucture_id,
                "resourceType": "BodyStructure",
                "identifier": [
                    {
                        "type": {
                            "coding": [
                                {
                                    "code": "tracking-uid",
                                    "display": "Tracking UID",
                                    "system": "http://hl7.org/fhir/uv/dicom-sr/CodeSystem/dicom-identifier-type",
                                }
                            ]
                        },
                        "use": "usual",
                        "value": structure_uid,
                    },
                    {
                        "type": {
                            "coding": [
                                {
                                    "code": "112040",
                                    "display": "Tracking Unique Identifier",
                                    "system": "http://dicom.nema.org/resources/ontology/DCM",
                                }
                            ]
                        },
                        "use": "usual",
                        "value": structure_uid,
                    },
                ],
                "location": {
                    "coding": [
                        {
                            "code": code.value,
                            "display": code.meaning,
                            "system": code_system,
                        }
                    ],
                },
                "morphology": {
                    "coding": [
                        {
                            "code": "91723000",
                            "display": "Anatomical Structure",
                            "system": "http://snomed.info/sct",
                        }
                    ],
                },
                "patient": {"reference": patient_ref},
            }
            if laterality is not None:
                json_data["locationQualifier"] = (
                    {
                        "coding": [
                            {
                                "code": code.value,
                                "display": code.meaning,
                                "system": code_system,
                            }
                        ],
                    },
                )

            element = construct_fhir_element("BodyStructure", json_data)
            body_structures.append(element)
        return body_structures

    @cached_property
    def _imaging_measurement_groups(self):
        """Create  FHIR ImagingMeasurements.

        Mapping the information from a
        DICOM TID 300 Measurement
        See: https://build.fhir.org/ig/HL7/dicom-sr/StructureDefinition-imaging-measurement.html

        Wrap each ImagingMeasurement into a FHIR ImagingMeasurementGroups,
        mapping the information from a DICOM TID 1411 Volumetric ROI Measurements and Qualitative Evaluations

        See: https://build.fhir.org/ig/HL7/dicom-sr/StructureDefinition-imaging-measurement-group.html

        """
        patient_ref = f"#{self._patient.id}"
        imaging_study_ref = f"#{self._imaging_study.id}"
        # device_ref = f"#{self._device.id}" # #NOTE: using the algorithm as the device, instead of the SR device name
        algorithm_identificaion_ref = f"#{self._algorithm_identification.id}"
        source_series_for_segmentation_ref = (
            f"#{self._source_series_for_segmentation.id}"
        )
        dcm_details = self._t1w_image_properties
        measurements = []
        for x, y, z in zip(
            self._ordered_seg_attrs,
            self._body_structures,
            self._referenced_segment_selections,
        ):
            (
                _,
                label,
                measure_code,
                _,
                structure_uid,
                measure_uid,
            ) = self._get_seg_attr_vals(x)
            measure_uid = x["measure_tracking_uid"]
            if measure_code.scheme_designator == "SCT":
                code_system = "http://snomed.info/sct"
            elif measure_code.scheme_designator == "DCM":
                code_system = "http://dicom.nema.org/resources/ontology/DCM"
            else:
                code_system = None  # #TODO: handle different systems
            body_structure_ref = f"#{y.id}"
            measurement_id = label.title().replace(" ", "").replace("_", "")
            vol = self.volume_measurements[x["segment_description"]]
            tid_300_measurement = {
                "id": f"{measurement_id}-VolumeMeasurement",
                "resourceType": "Observation",
                "identifier": [
                    {
                        "type": {
                            "coding": [
                                {
                                    "code": "observation-uid",
                                    "display": "Observation UID",
                                    "system": "http://hl7.org/fhir/uv/dicom-sr/CodeSystem/dicom-identifier-type",
                                }
                            ]
                        },
                        "use": "usual",
                        "value": measure_uid,
                    },
                    {
                        "type": {
                            "coding": [
                                {
                                    "code": "112040",
                                    "display": "Tracking Unique Identifier",
                                    "system": "http://dicom.nema.org/resources/ontology/DCM",
                                }
                            ]
                        },
                        "use": "usual",
                        "value": measure_uid,
                    },
                ],
                "partOf": [{"reference": imaging_study_ref}],
                "category": [
                    {
                        "coding": [
                            {
                                "code": "125007",
                                "display": "Measurement Group",
                                "system": "http://dicom.nema.org/resources/ontology/DCM",
                            }
                        ],
                    },
                ],
                "code": {
                    "coding": [
                        {
                            "code": "68604-8",
                            "display": "Radiology Diagnostic study note",
                            "system": "http://loinc.org",
                        }
                    ]
                },
                "subject": {"reference": patient_ref},
                # #TODO: should both bodySite and focus be used for same coding?
                "focus": [{"reference": body_structure_ref}],
                "issued": observation_datetime,
                "valueQuantity": {
                    "value": vol,
                    "unit": "mm3",
                    "system": "http://unitsofmeasure.org/",
                },
                # NOTE: The spec uses bodyStructure, not yet implemented (Reference(BodyStructure) resource)
                # https://build.fhir.org/ig/HL7/dicom-sr/StructureDefinition-imaging-measurement.html
                "bodySite": {
                    "coding": [
                        {
                            "code": measure_code.value,
                            "display": measure_code.meaning,
                            "system": code_system,
                        }
                    ]
                },
                "device": {"reference": algorithm_identificaion_ref},
                "effectiveDateTime": dcm_details["study_datetime"],
                "status": "final",
                # #TODO: referenced range can be used for norm stats
                # "referenced_range": [...]
            }

            tid_300_element = construct_fhir_element("Observation", tid_300_measurement)

            referenced_segment_ref = f"#{z.id}"
            tid_300_element_ref = f"#{tid_300_element.id}"
            tid_1411_measurement_group = {
                "id": f"{measurement_id}-VolumeMeasurementGroup",
                "resourceType": "Observation",
                "identifier": [
                    {
                        "type": {
                            "coding": [
                                {
                                    "code": "observation-uid",
                                    "display": "Observation UID",
                                    "system": "http://hl7.org/fhir/uv/dicom-sr/CodeSystem/dicom-identifier-type",
                                }
                            ]
                        },
                        "use": "usual",
                        "value": structure_uid,
                    },
                    {
                        "type": {
                            "coding": [
                                {
                                    "code": "112040",
                                    "display": "Tracking Unique Identifier",
                                    "system": "http://dicom.nema.org/resources/ontology/DCM",
                                }
                            ]
                        },
                        "use": "usual",
                        "value": structure_uid,
                    },
                ],
                "partOf": tid_300_measurement["partOf"],
                "category": tid_300_measurement["category"],
                "code": tid_300_measurement["code"],
                "subject": tid_300_measurement["subject"],
                "focus": tid_300_measurement["focus"]
                + [
                    {"reference": referenced_segment_ref},
                    {"reference": source_series_for_segmentation_ref},
                ],
                "contained": [tid_300_element],
                "status": "final",
                "hasMember": [{"reference": tid_300_element_ref}],
            }
            tid_1411_element = construct_fhir_element(
                "Observation", tid_1411_measurement_group
            )
            # #TODO figure out what to return, or whether to contains?
            measurements.append(tid_1411_element)
        return measurements

    # #NOTE: no longer using until ImagingSelection is implemented
    # @cached_property
    # def _referenced_segment_instances(self):
    #     referenced_seg_resources = []
    #     for x, y in zip(self._ordered_seg_attrs, self._referenced_segments):
    #         (num, label, _, _, _, _) = self._get_seg_attr_vals(x)
    #         referenced_sop_sequence = y.ReferencedSOPSequence[0]
    #         assert num == referenced_sop_sequence.ReferencedSegmentNumber
    #         instance_uid = referenced_sop_sequence.ReferencedSOPInstanceUID
    #         json_data = {
    #             "id": f"Instance-ReferencedSegment-{label.title().replace(' ', '').replace('_', '')}",
    #             "resourceType": "Observation",
    #             "uid": instance_uid,
    #             "sopClass": {
    #                 "code": "urn:oid:1.2.840.10008.5.1.4.1.1.66.4",
    #                 "display": "Segmentation Storage",
    #                 "system": "urn:ietf:rfc:3986",
    #             },
    #             "subset": [str(num)],
    #         }
    #         element = construct_fhir_element("Observation", json_data)
    #         referenced_seg_resources.append(element)
    #     return referenced_seg_resources

    @cached_property
    def _referenced_segment_selections(self):
        imaging_study_ref = f"#{self._imaging_study.id}"
        # patient_ref = f"#{self._patient.id}"
        imaging_study_ref = f"#{self._imaging_study.id}"
        referenced_seg_selections = []

        for x, y in zip(self._ordered_seg_attrs, self._referenced_segments):
            (num, label, _, _, _, _) = self._get_seg_attr_vals(x)
            referenced_sop_sequence = y.ReferencedSOPSequence[0]
            assert num == referenced_sop_sequence.ReferencedSegmentNumber
            instance_uid = referenced_sop_sequence.ReferencedSOPInstanceUID
            json_data = {
                "id": f"ImagingSelection-ReferencedSegment-{label.title().replace(' ', '').replace('_', '')}",
                "resourceType": "Observation",
                "code": {
                    "coding": [
                        {
                            "code": "121191",
                            "display": "Referenced Segment",
                            "system": "http://dicom.nema.org/resources/ontology/DCM",
                        }
                    ]
                },
                "derivedFrom": [{"reference": imaging_study_ref}],
                "status": "available",
                # "instance": [orjson.loads(y.json())],
                "component": [
                    {
                        "valueCodeableConcept": {
                            "coding": [
                                {
                                    "code": "urn:oid:1.2.840.10008.5.1.4.1.1.66.4",
                                    "display": "Segmentation Storage",
                                    "system": "urn:ietf:rfc:3986",
                                }
                            ]
                        },
                        "code": {
                            "coding": [
                                {
                                    "code": "110181",
                                    "display": "SOP Class UID",
                                    "system": "http://dicom.nema.org/resources/ontology/DCM",
                                }
                            ]
                        },
                    },
                    {
                        "valueString": instance_uid,
                        "code": {
                            "coding": [
                                {
                                    "code": "referenced-sop-instance-uid",
                                    "display": "Referenced SOP Instance UID",
                                    "system": "http://hl7.org/fhir/uv/dicom-sr/CodeSystem/dicom-identifier-type",
                                }
                            ]
                        },
                    },
                ],
                "valueString": str(num),
            }
            element = construct_fhir_element("Observation", json_data)
            referenced_seg_selections.append(element)
        return referenced_seg_selections

    @cached_property
    def _source_series_for_segmentation(self):
        patient_ref = f"#{self._patient.id}"
        imaging_study_ref = f"#{self._imaging_study.id}"
        dcm_details = self._t1w_image_properties

        json_data = {
            "id": "ImagingSelection-SourceSeriesForSegmentation",
            "resourceType": "Observation",
            "subject": {"reference": patient_ref},
            "code": {
                "coding": [
                    {
                        "code": "121232",
                        "display": "Source series for segmentation",
                        "system": "http://dicom.nema.org/resources/ontology/DCM",
                    }
                ]
            },
            "derivedFrom": [{"reference": imaging_study_ref}],
            # "seriesUid": dcm_details["SeriesInstanceUID"],
            "valueCodeableConcept": {
                "coding": [
                    {
                        "code": "Series Instance UID",
                        "display": "Series Instance UID",
                        "system": "http://dicom.nema.org/resources/ontology/DCM",
                    },
                ],
                "text": dcm_details["SeriesInstanceUID"],
            },
            "status": "available",
        }
        element = construct_fhir_element("Observation", json_data)
        return element

    @cached_property
    def _tid1411_vol_rois(self):
        measurements = hd.sr.utils.find_content_items(
            dataset=self.sr, name=codes.DCM.ImagingMeasurements, recursive=True
        )
        return measurements[0]

    @cached_property
    def _referenced_segments(self):
        return hd.sr.utils.find_content_items(
            dataset=self._tid1411_vol_rois,
            name=codes.DCM.ReferencedSegment,
            recursive=True,
        )

    @cached_property
    def _t1w_image_properties(self):
        dcm = self.t1w_image_dataset
        dcm_attrs = [
            x
            for x in dir(dcm)
            if not x.startswith("_") and not inspect.ismethod(getattr(dcm, x))
        ]
        dcm_properties = {x: getattr(dcm, x) for x in dcm_attrs}

        # add in extra keywords
        dcm_properties[
            "study_datetime"
        ] = f"{pydicom.valuerep.DA(dcm.StudyDate).isoformat()}T{pydicom.valuerep.TM(dcm.StudyTime).isoformat()}"
        dcm_properties["num_instances"] = len(self.t1w_image_datasets)
        # research scans may not have AccessionNumbers
        try:
            acsn = dcm.AccessionNumber
        except AttributeError:
            acsn = "No AccessionNumber provided"
        if acsn == "":
            acsn = "No AccessionNumber provided"
        dcm_properties["acsn"] = acsn
        return dcm_properties

    @cached_property
    def fhir(self):
        report_id = generate_uid(CMET_ROOT_UID)
        dcm_details = self._t1w_image_properties

        imaging_study_ref = f"#{self._imaging_study.id}"

        imaging_measurement_groups_refs = [
            {"reference": f"#{x.id}"} for x in self._imaging_measurement_groups
        ]

        dcm_details = self._t1w_image_properties
        json_data = {
            "id": report_id,
            "resourceType": "DiagnosticReport",
            "code": {
                "coding": [
                    {
                        "code": "68604-8",
                        "display": "Radiology Diagnostic study note",
                        "system": "http://loinc.org",
                    }
                ]
            },
            "category": [
                {
                    "coding": [
                        {
                            "code": "RAD",
                            "display": "Radiology",
                            "system": "http://terminology.hl7.org/CodeSystem/v2-0074",
                        }
                    ]
                }
            ],
            "effectiveDateTime": dcm_details["study_datetime"],
            "issued": observation_datetime,
            "result": imaging_measurement_groups_refs,
            "imagingStudy": [{"reference": imaging_study_ref}],
            "status": "final",
            "conclusion": REPORT_NAME,
            "text": {
                "status": "generated",
                "div": f'<div xmlns="http://www.w3.org/1999/xhtml"><p><b>{REPORT_NAME}</b></p></div>',
            },
            "contained": [
                self._patient,
                self._imaging_study,
                self._device,
                self._algorithm_identification,
                self._source_series_for_segmentation,
                *self._referenced_segment_selections,
                *self._body_structures,
                *self._imaging_measurement_groups,
            ],
        }
        element = construct_fhir_element("DiagnosticReport", json_data)

        if self.fhir_output_file is not None:
            with open(self.fhir_output_file, "w", encoding="utf8") as f:
                f.write(element.json(indent=2))

        return element

    def to_df(self, output_filename: Union[str, Path, None] = None) -> pd.DataFrame:
        df = pd.DataFrame(
            self.volume_measurements.items(), columns=["Structure", "Volume"]
        )
        # add in other useful information
        for x in self._seg_attrs:
            struct = x["segment_description"]
            for y in [
                "label_id",
                "recommended_rgb_value",
                "segment_description",
                "segmented_property_type",
                "segmented_property_category",
                "segmented_property_modifier",
            ]:
                try:
                    df.loc[df["Structure"] == struct, y] = str(x[y])
                except KeyError:
                    pass
        if output_filename is not None:
            df.to_csv(output_filename, index=False)
        return df


# pylint: disable=too-many-lines
