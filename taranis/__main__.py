import logging

import click
import concurrent.futures
import glob
import rich.console
import rich.logging
import rich.traceback
import sys
import time

import taranis.utils
import taranis.analyze_schema
import taranis.reference_alleles
import taranis.allele_calling

import taranis.inferred_alleles

log = logging.getLogger()

# Set up rich stderr console
stderr = rich.console.Console(
    stderr=True, force_terminal=taranis.utils.rich_force_colors()
)


def expand_wildcards(ctx, param, value):
    if value:
        expanded_paths = []
        for path in value:
            # Check if path contains wildcard
            if "*" in path:
                # Expand wildcard
                expanded_paths.extend(glob.glob(path))
            else:
                expanded_paths.append(path)
        return expanded_paths
    return None


def run_taranis():
    # Set up the rich traceback
    rich.traceback.install(console=stderr, width=200, word_wrap=True, extra_lines=1)

    # Print taranis header
    stderr.print(
        "[blue]                ______           ___                     ___    ",
        highlight=False,
    )
    stderr.print(
        "[blue]   \    |-[grey39]-|  [blue] |__--__|   /\    |   \    /\    |\   | | |  ",
        highlight=False,
    )
    stderr.print(
        "[blue]    \   \  [grey39]/ [blue]     ||     /  \   |__ /   /  \   | \  | | |___   ",
        highlight=False,
    )
    stderr.print(
        "[blue]    /  [grey39] / [blue] \      ||    /____\  |  \   /____\  |  \ | |     |",
        highlight=False,
    )
    stderr.print(
        "[blue]   /   [grey39] |-[blue]-|      ||   /      \ |   \ /      \ |   \| |  ___| ",
        highlight=False,
    )

    # stderr.print("[green]                                          `._,._,'\n", highlight=False)
    __version__ = "3.0.0"
    stderr.print(
        "\n" "[grey39]    Taranis version {}".format(__version__), highlight=False
    )
    # Lanch the click cli
    taranis_cli()


# Customise the order of subcommands for --help
class CustomHelpOrder(click.Group):
    def __init__(self, *args, **kwargs):
        self.help_priorities = {}
        super(CustomHelpOrder, self).__init__(*args, **kwargs)

    def get_help(self, ctx):
        self.list_commands = self.list_commands_for_help
        return super(CustomHelpOrder, self).get_help(ctx)

    def list_commands_for_help(self, ctx):
        """reorder the list of commands when listing the help"""
        commands = super(CustomHelpOrder, self).list_commands(ctx)
        return (
            c[1]
            for c in sorted(
                (self.help_priorities.get(command, 1000), command)
                for command in commands
            )
        )

    def command(self, *args, **kwargs):
        """Behaves the same as `click.Group.command()` except capture
        a priority for listing command names in help.
        """
        help_priority = kwargs.pop("help_priority", 1000)
        help_priorities = self.help_priorities

        def decorator(f):
            cmd = super(CustomHelpOrder, self).command(*args, **kwargs)(f)
            help_priorities[cmd.name] = help_priority
            return cmd

        return decorator


@click.group(cls=CustomHelpOrder)
@click.version_option(taranis.__version__)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    default=False,
    help="Print verbose output to the console.",
)
@click.option(
    "-l", "--log-file", help="Save a verbose log to a file.", metavar="filename"
)
def taranis_cli(verbose, log_file):
    # Set the base logger to output DEBUG
    log.setLevel(logging.DEBUG)

    # Set up logs to a file if we asked for one
    if log_file:
        log_fh = logging.FileHandler(log_file, encoding="utf-8")
        log_fh.setLevel(logging.DEBUG)
        log_fh.setFormatter(
            logging.Formatter(
                "[%(asctime)s] %(name)-20s [%(levelname)-7s]  %(message)s"
            )
        )
        log.addHandler(log_fh)


@taranis_cli.command(help_priority=1)
@click.option(
    "-i",
    "--inputdir",
    required=True,
    multiple=False,
    type=click.Path(),
    help="Directory where the schema with the core gene files are located. ",
)
@click.option(
    "-o",
    "--output",
    required=True,
    multiple=False,
    type=click.Path(),
    help="Output folder to save analyze schema",
)
@click.option(
    "--remove-subset/--no-remove-subset",
    required=False,
    default=False,
    help="Remove allele subsequences from the schema.",
)
@click.option(
    "--remove-duplicated/--no-remove-duplicated",
    required=False,
    default=False,
    help="Remove duplicated subsequences from the schema.",
)
@click.option(
    "--remove-no-cds/--no-remove-no-cds",
    required=False,
    default=False,
    help="Remove no CDS alleles from the schema.",
)
@click.option(
    "--output-allele-annot/--no-output-allele-annot",
    required=False,
    default=True,
    help="output prokka/allele annotation for all alleles in locus",
)
@click.option(
    "--genus",
    required=False,
    default="Genus",
    help="Genus name for Prokka schema genes annotation. Default is Genus.",
)
@click.option(
    "--species",
    required=False,
    default="species",
    help="Species name for Prokka schema genes annotation. Default is species",
)
@click.option(
    "--usegenus",
    required=False,
    default="Genus",
    help="Use genus-specific BLAST databases for Prokka schema genes annotation (needs --genus). Default is False.",
)
@click.option(
    "--cpus",
    required=False,
    multiple=False,
    type=int,
    default=1,
    help="Number of cpus used for execution",
)
def analyze_schema(
    inputdir: str,
    output: str,
    remove_subset: bool,
    remove_duplicated: bool,
    remove_no_cds: bool,
    output_allele_annot: bool,
    genus: str,
    species: str,
    usegenus: str,
    cpus: int,
):
    _ = taranis.utils.check_additional_programs_installed([["prokka", "--version"]])
    schema_files = taranis.utils.get_files_in_folder(inputdir, "fasta")

    results = []
    max_cpus = taranis.utils.cpus_available()
    if cpus > max_cpus:
        stderr.print("[red] Number of CPUs bigger than the CPUs available")
        stderr.print("Running code with ", max_cpus)
        cpus = max_cpus
    # Keeping 3 threads for running prokka for each parallel process
    using_cpus, prokka_cpus = [cpus // 3, 3] if cpus // 3 >= 1 else [1, 1]
    start = time.perf_counter()
    with concurrent.futures.ThreadPoolExecutor(max_workers=using_cpus) as executor:
        futures = [
            executor.submit(
                taranis.analyze_schema.parallel_execution,
                schema_file,
                output,
                remove_subset,
                remove_duplicated,
                remove_no_cds,
                genus,
                species,
                usegenus,
                prokka_cpus,
            )
            for schema_file in schema_files
        ]
        # Collect results as they complete
        for future in concurrent.futures.as_completed(futures):
            results.append(future.result())
    _ = taranis.analyze_schema.collect_statistics(results, output, output_allele_annot)

    finish = time.perf_counter()
    print(f"Schema analyze finish in {round((finish-start)/60, 2)} minutes")


# Reference alleles
@taranis_cli.command(help_priority=2)
@click.option(
    "-s",
    "--schema",
    required=True,
    multiple=False,
    type=click.Path(),
    help="Directory where the schema with the core gene files are located. ",
)
@click.option(
    "-o",
    "--output",
    required=True,
    multiple=False,
    type=click.Path(),
    help="Output folder to save reference alleles",
)
@click.option(
    "--eval-cluster/--no-eval-cluster",
    required=False,
    default=True,
    help="Evaluate if the reference alleles match against blast with a 90% identity",
)
@click.option(
    "-k",
    "--kmer-size",
    required=False,
    type=int,
    default=21,
    help="Mash parameter for K-mer size.",
)
@click.option(
    "-S",
    "--sketch-size",
    required=False,
    type=int,
    default=2000,
    help="Mash parameter for Sketch size",
)
@click.option(
    "-r",
    "--cluster-resolution",
    required=False,
    type=float,
    default=0.75,
    help="Resolution value used for clustering.",
)
@click.option(
    "--seed",
    required=False,
    type=int,
    default=None,
    help="Seed value for clustering",
)
@click.option(
    "--cpus",
    required=False,
    multiple=False,
    type=int,
    default=1,
    help="Number of cpus used for execution",
)
@click.option(
    "--force/--no-force",
    required=False,
    default=False,
    help="Overwrite the output folder if it exists",
)
def reference_alleles(
    schema: str,
    output: str,
    eval_cluster: bool,
    kmer_size: int,
    sketch_size: int,
    cluster_resolution: float,
    seed: int,
    cpus: int,
    force: bool,
):
    _ = taranis.utils.check_additional_programs_installed(
        [["mash", "--version"], ["makeblastdb", "-version"], ["blastn", "-version"]]
    )
    start = time.perf_counter()
    max_cpus = taranis.utils.cpus_available()
    if cpus > max_cpus:
        stderr.print("[red] Number of CPUs bigger than the CPUs available")
        stderr.print("Running code with ", max_cpus)
        cpus = max_cpus
    schema_files = taranis.utils.get_files_in_folder(schema, "fasta")

    # Check if output folder exists
    if not force:
        _ = taranis.utils.prompt_user_if_folder_exists(output)
    """Create the reference alleles from the schema """
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=cpus) as executor:
        futures = [
            executor.submit(
                taranis.reference_alleles.parallel_execution,
                f_file,
                output,
                eval_cluster,
                kmer_size,
                sketch_size,
                cluster_resolution,
                seed,
            )
            for f_file in schema_files
        ]
        for future in concurrent.futures.as_completed(futures):
            try:
                results.append(future.result())
            except Exception as e:
                print(e)
                continue
    _ = taranis.reference_alleles.collect_statistics(results, eval_cluster, output)
    finish = time.perf_counter()
    print(f"Reference alleles finish in {round((finish-start)/60, 2)} minutes")


@taranis_cli.command(help_priority=3)
@click.option(
    "-s",
    "--schema",
    required=True,
    multiple=False,
    type=click.Path(exists=True),
    help="Directory where the schema with the core gene files are located. ",
)
@click.option(
    "-r",
    "--reference",
    required=True,
    multiple=False,
    type=click.Path(exists=True),
    help="Directory where the schema reference allele files are located. ",
)
@click.option(
    "-a",
    "--annotation",
    required=True,
    multiple=False,
    type=click.Path(exists=True),
    help="Annotation file. ",
)
@click.option(
    "-t",
    "--threshold",
    required=False,
    nargs=1,
    default=0.8,
    type=float,
    help="Threshold value to consider in blast. Values from 0 to 1. default 0.8",
)
@click.option(
    "-p",
    "--perc-identity",
    required=False,
    nargs=1,
    default=90,
    type=int,
    help="Percentage of identity to consider in blast. default 90",
)
@click.option(
    "-o",
    "--output",
    required=True,
    multiple=False,
    type=click.Path(),
    help="Output folder to save reference alleles",
)
@click.option(
    "--force/--no-force",
    required=False,
    default=False,
    help="Overwrite the output folder if it exists",
)
@click.argument(
    "assemblies",
    callback=expand_wildcards,
    nargs=-1,
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    "--snp/--no-snp",
    required=False,
    default=False,
    help="Create SNP file for alleles in assembly in relation with reference allele",
)
@click.option(
    "--alignment/--no-alignment",
    required=False,
    default=False,
    help="Create alignment files",
)
@click.option(
    "-q",
    "--proteine-threshold",
    required=False,
    nargs=1,
    default=80,
    type=int,
    help="Threshold of protein coverage to consider as TPR. default 90",
)
@click.option(
    "-i",
    "--increase-sequence",
    required=False,
    nargs=1,
    default=20,
    type=int,
    help="Increase the number of triplet sequences to find the stop codon. default 20",
)
@click.option(
    "--cpus",
    required=False,
    multiple=False,
    type=int,
    default=1,
    help="Number of cpus used for execution",
)
def allele_calling(
    schema: str,
    reference: str,
    annotation: str,
    assemblies: list,
    threshold: float,
    perc_identity: int,
    output: str,
    force: bool,
    snp: bool,
    alignment: bool,
    proteine_threshold: int,
    increase_sequence: int,
    cpus: int,
):
    _ = taranis.utils.check_additional_programs_installed(
        [["blastn", "-version"], ["makeblastdb", "-version"], ["mafft", "--version"]]
    )
    schema_ref_files = taranis.utils.get_files_in_folder(reference, "fasta")
    if len(schema_ref_files) == 0:
        log.error("Referenc allele folder %s does not have any fasta file", schema)
        stderr.print("[red] reference allele folder does not have any fasta file")
        sys.exit(1)

    # Check if output folder exists
    if not force:
        _ = taranis.utils.prompt_user_if_folder_exists(output)
    # Filter fasta files from reference folder
    # ref_alleles = glob.glob(os.path.join(reference, "*.fasta"))
    max_cpus = taranis.utils.cpus_available()
    if cpus > max_cpus:
        stderr.print("[red] Number of CPUs bigger than the CPUs available")
        stderr.print("Running code with ", max_cpus)
        cpus = max_cpus
    # Read the annotation file
    stderr.print("[green] Reading annotation file")
    log.info("Reading annotation file")
    map_pred = [["gene", 7], ["product", 8], ["allele_quality", 9]]
    prediction_data = taranis.utils.read_compressed_file(
        annotation, separator=",", index_key=1, mapping=map_pred
    )
    # Create the instanace for inference alleles
    inf_allele_obj = taranis.inferred_alleles.InferredAllele()
    """Analyze the sample file against schema to identify alleles
    """

    start = time.perf_counter()
    results = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=cpus) as executor:
        futures = [
            executor.submit(
                taranis.allele_calling.parallel_execution,
                assembly_file,
                schema,
                prediction_data,
                schema_ref_files,
                threshold,
                perc_identity,
                output,
                inf_allele_obj,
                snp,
                alignment,
                proteine_threshold,
                increase_sequence,
            )
            for assembly_file in assemblies
        ]
        for future in concurrent.futures.as_completed(futures):
            try:
                results.append(future.result())
            except Exception as e:
                print(e)
                continue
    """
    for assembly_file in assemblies:
        results.append(
            taranis.allele_calling.parallel_execution(
                assembly_file,
                schema,
                prediction_data,
                schema_ref_files,
                threshold,
                perc_identity,
                output,
                inf_allele_obj,
                snp,
                alignment,
                proteine_threshold,
                increase_sequence,
            )
        )
    """
    _ = taranis.allele_calling.collect_data(
        results, output, snp, alignment, schema_ref_files, cpus
    )
    finish = time.perf_counter()
    print(f"Allele calling finish in {round((finish-start)/60, 2)} minutes")
    log.info("Allele calling finish in %s minutes", round((finish - start) / 60, 2))
    # sample_allele_obj.analyze_sample()
