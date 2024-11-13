use anyhow::Result;
use clap::Parser;
use log::{debug, info, LevelFilter, trace};

mod cli;
mod ava;
mod twoset;

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

    let estimates: Vec<(&[u8], f32)> = if let Some(num) = args.num_reads {
        info!("Running with all-vs-all strategy with {} reads", num);
        ava_strategy(args.input, num)?
    } else if let (Some(target_num_reads), Some(query_num_reads)) =
        (args.target_num_reads, args.query_num_reads)
    {
        twoset_strategy(args.input, target_num_reads, query_num_reads)?
    } else {
        unreachable!("No strategy could be determined. Please raise an issue at <https://github.com/mbhall88/lrge/issues>")
    };

    if log::log_enabled!(log::Level::Trace) {
        for (rid, est) in estimates {
            trace!("Estimate for {}: {}", String::from_utf8_lossy(rid), est);
        }
    }

    todo!("Determine the median of the estimates, depending on whether infinite estimates are to be included");
    let estimate = 0;
    println!("{}", estimate);
    Ok(())
}
