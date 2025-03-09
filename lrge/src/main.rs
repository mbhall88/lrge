use crate::utils::{create_temp_dir, format_estimate};
use anyhow::{bail, Context, Result};
use clap::Parser;
use liblrge::Estimate;
use log::{debug, info, LevelFilter};
use std::fs::File;
use std::io;
use std::io::Write;

mod cli;
mod utils;

fn setup_logging(quiet: u8, verbose: u8) {
    let sum = verbose as i8 - quiet as i8;

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
        info!(
            "Created temporary directory at {}",
            tmpdir.path().to_string_lossy()
        );
    } else {
        debug!(
            "Created temporary directory at {}",
            tmpdir.path().to_string_lossy()
        );
    }

    let mut output: Box<dyn Write> = if args.output == "-" {
        Box::new(io::stdout())
    } else {
        Box::new(File::create(&args.output).context("Failed to create output file")?)
    };

    let mut strategy: Box<dyn Estimate> = if let Some(num) = args.num_reads {
        info!("Running all-vs-all strategy with {} reads", num);
        let builder = liblrge::ava::Builder::new()
            .num_reads(num)
            .remove_internal(args.do_filt, args.max_overhang_ratio)
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
            .remove_internal(args.do_filt, args.max_overhang_ratio)
            .threads(args.threads)
            .tmpdir(tmpdir.path())
            .seed(args.seed);

        Box::new(builder.build(args.input))
    } else {
        unreachable!("No strategy could be determined. Please raise an issue at <https://github.com/mbhall88/lrge/issues>")
    };

    let est_result = strategy
        .estimate(!args.with_infinity, Some(args.lower_q), Some(args.upper_q))
        .context("Failed to generate estimate")?;

    let estimate = est_result.estimate;
    let low_q = est_result.lower;
    let upper_q = est_result.upper;

    match estimate {
        Some(est) => {
            let formatted_est = format_estimate(est);
            let mut msg = format!("Estimated genome size: {formatted_est}");
            if let (Some(low), Some(high)) = (low_q, upper_q) {
                let formatted_low = format_estimate(low);
                let formatted_high = format_estimate(high);
                msg.push_str(&format!(" (IQR: {formatted_low} - {formatted_high})"));
            }
            info!("{}", msg);

            if args.precise {
                writeln!(output, "{est}")?;
            } else {
                writeln!(output, "{est:.0}")?;
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
