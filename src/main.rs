use clap::Parser;
use log::{debug, LevelFilter};

mod cli;

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

fn main() {
    let opts = cli::Opts::parse();
    setup_logging(opts.quiet, opts.verbose);
    debug!("{:?}", opts)
}
