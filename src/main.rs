use std::env;
use std::io;
use std::io::Write;
use std::process::exit;

mod taylormaccoll;
mod busemann;
mod inlet;
mod utils;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("usage: {} [inlet type]", args[0]);
        std::process::exit(1);
    }

    let inlet_type = &args[1];

    match inlet_type.as_str() {
        "b" | "busemann" => {
            // ask use to input inlet mach number
            print!("enter the design exit mach number: ");
            io::stdout().flush().unwrap();
            
            let mut input: String = String::new();
            io::stdin().read_line(&mut input)
                .expect("failed to read input mach");
            let exit_mach: f64 = match input.trim().parse() {
                Ok(num) => num,
                Err(_) => {
                    eprintln!("invalid exit mach number");
                    exit(1);
                }
            };
            print!("enter the design free stream mach number: ");
            io::stdout().flush().unwrap();

            input.clear();

            io::stdin().read_line(&mut input)
                .expect("failed to read input");
            let freestream_mach: f64 = match input.trim().parse() {
                Ok(num) => num,
                Err(_) => {
                    eprintln!("invalid freestream mach number");
                    exit(1);
                }
            };
            println!("{}, {}", exit_mach, freestream_mach);
            let busemann: inlet::Inlet = busemann::calc_contour(exit_mach, freestream_mach);
            busemann.export_csv();
            busemann.plot("busemann.png");
        }
        "icfa" => {
            todo!("icfa inlet")
        }
        "tb" | "truncated busemann" => {
            todo!("truncated busemann inlet")
        }
        _ => {
            eprintln!("unknown inlet type '{}'", inlet_type);
            exit(1);
        }
    }
}
