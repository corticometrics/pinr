# type: ignore[attr-defined]
from typing import Optional

import typer
from rich.console import Console

from pinr import __version__
from pinr.model import FreeSurferSeg, template_file

app = typer.Typer(
    name="pinr",
    help="Python tools for Interoperable Neuromorphometry Reporting",
    add_completion=True,
)
console = Console()


def version_callback(print_version: bool) -> None:
    """Print the version of the package."""
    if print_version:
        console.print(f"[yellow]pinr[/] version: [bold blue]{__version__}[/]")
        raise typer.Exit()


@app.command()
def main(
    t1w_dicom_file: str = typer.Option(
        ...,
        "-i",
        "--t1w",
        help="One T1w DICOM file, part of the series used as input to FreeSurfer.",
        show_default=False,
    ),
    aseg_file: str = typer.Option(
        ...,
        "-s",
        "--aseg",
        "--seg",
        help="aseg.mgz Segmentation file created by FreeSurfer",
        show_default=False,
    ),
    aseg_stats_file: str = typer.Option(
        ...,
        "-m",
        "--measures",
        "--stats",
        help="aseg.stats file output by FreeSurfer.",
        show_default=False,
    ),
    dicom_seg_output_file: Optional[str] = typer.Option(
        None, "--dcm-seg", help="output file name for DICOM Seg", show_default=False
    ),
    dicom_sr_output_file: Optional[str] = typer.Option(
        None, "--dcm-sr", help="output file name for DICOM SR.", show_default=False
    ),
    fhir_output_file: Optional[str] = typer.Option(
        None,
        "--fhir",
        help="output file name for FHIR Diagnostic Report JSON",
        show_default=False,
    ),
    dicom_seg_schema_file: Optional[str] = typer.Option(
        str(template_file),
        "--schema",
        help="Schema file used to provided coded inputs for each structure, by default uses the one provided",
    ),
    print_version: bool = typer.Option(
        None,
        "-v",
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the pinr package.",
    ),
) -> None:
    # pylint: disable=unused-argument
    """Convert FreeSurfer aseg results into medical imaging formats."""
    console.print("[bold]Initializing FreeSurfer Seg conversion model...[/]")
    aseg = FreeSurferSeg(
        t1w_dicom_file=t1w_dicom_file,
        aseg_file=aseg_file,
        aseg_stats_file=aseg_stats_file,
        dicom_seg_output_file=dicom_seg_output_file,
        dicom_sr_output_file=dicom_sr_output_file,
        fhir_output_file=fhir_output_file,
        dicom_seg_schema_file=dicom_seg_schema_file,
    )
    if dicom_seg_output_file is not None:
        console.print(f"[bold]Saving DICOM Seg to: {dicom_seg_output_file}[/]")
        _ = aseg.seg
    if dicom_sr_output_file is not None:
        console.print(f"[bold]Saving DICOM SR to: {dicom_sr_output_file}[/]")
        _ = aseg.sr
    if fhir_output_file is not None:
        console.print(
            f"[bold]Saving FHIR Diagnostic Report JSON to: {fhir_output_file}[/]"
        )
        _ = aseg.fhir


if __name__ == "__main__":
    app()
