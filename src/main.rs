use clap::Parser;

mod cli;

fn main() {
    let opts = cli::Opts::parse();
    println!("{:?}", opts)
}
