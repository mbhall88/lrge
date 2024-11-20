use crate::utils::{create_temp_dir, format_estimate};
use anyhow::{bail, Context, Result};
use clap::Parser;
use liblrge::Estimate;
use log::{debug, info, LevelFilter};

mod cli;
mod utils;

fn setup_logging(quiet: u8, verbose: u8) {
    let sum = (verbose - quiet) as i16;
    let lvl = match sum {
        1 => LevelFilter::Debug,
        2.. => LevelFilter::Trace,
        -1 => LevelFilter::Warn,
        -2 => LevelFilter::Error,
        i if i < -2 => LevelFilter::Off,
        _ => LevelFilter::Info,
    };
    let mut log_builder = env_logger::Builder::new();
    log_builder
        .filter(None, lvl)
        .filter_module("mio", LevelFilter::Off)
        .filter_module("reqwest", LevelFilter::Off);
    log_builder.init();
}

fn main() -> Result<()> {
    let args = cli::Args::parse();
    setup_logging(args.quiet, args.verbose);
    debug!("{:?}", args);

    let tmpdir = create_temp_dir(args.temp_dir.as_ref(), args.keep_temp)?;
    if args.keep_temp {
        info!("Created temporary directory at {:?}", tmpdir.path());
    } else {
        debug!("Created temporary directory at {:?}", tmpdir.path());
    }

    let mut strategy: Box<dyn Estimate> = if let Some(num) = args.num_reads {
        info!("Running all-vs-all strategy with {} reads", num);
        let builder = liblrge::ava::Builder::new()
            .num_reads(num)
            .threads(args.threads)
            .tmpdir(tmpdir.path())
            .seed(args.seed);

        Box::new(builder.build(args.input))
    } else if let (Some(target_num_reads), Some(query_num_reads)) =
        (args.target_num_reads, args.query_num_reads)
    {
        info!(
            "Running two-set strategy with {} target reads and {} query reads",
            target_num_reads, query_num_reads
        );
        let builder = liblrge::twoset::Builder::new()
            .target_num_reads(target_num_reads)
            .query_num_reads(query_num_reads)
            .threads(args.threads)
            .tmpdir(tmpdir.path())
            .seed(args.seed);

        Box::new(builder.build(args.input))
    } else {
        unreachable!("No strategy could be determined. Please raise an issue at <https://github.com/mbhall88/lrge/issues>")
    };

    let estimate = if args.with_infinity {
        strategy.estimate_with_infinity()
    } else {
        strategy.estimate()
    }
    .context("Failed to generate estimate")?;

    match estimate {
        Some(est) => {
            let formatted_est = format_estimate(est);
            info!("Estimated genome size: {formatted_est}");

            if args.precise {
                println!("{est}");
            } else {
                println!("{est:.0}");
            }
        }
        None => {
            if args.with_infinity {
                bail!("No estimates were generated")
            } else {
                bail!("No finite estimates were generated")
            }
        }
    }

    info!("Done!");
    Ok(())
}
