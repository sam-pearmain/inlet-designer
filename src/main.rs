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

    let inlet_type: &String = &args[1];

    match inlet_type.as_str() {
        "b" | "busemann" => {
            print!(
                "please select design methodology:\n
                 - moc freestream & exit mach [1]\n
                 - moc exit mach & compression efficiency [2]\n
                 - taylor maccoll freestream & exit mach [3]\n
                 - taylor maccoll exit mach & compression efficiency [4]\n: "
            );
            io::stdout().flush().unwrap();

            let mut input = String::new();
            io::stdin().read_line(&mut input)
                .expect("failed to read input method");
            match input.trim().parse().unwrap() {
                1 => todo!(),
                2 => todo!(),
                3 => {
                    // ask use to input inlet mach number
                    print!("enter the design exit mach number: ");
                    io::stdout().flush().unwrap();

                    input.clear();
                    
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
                4 => todo!(),
                _ => panic!("unknown method for designing busemann inlet, select [1], [2], [3], or [4]")
            }
        }
        "icfa" => {
            todo!("icfa inlet")
        }
        "tb" | "truncated-busemann" => {
            todo!("truncated busemann inlet")
        }
        "bcb" | "boundary-corrected-busemann" => {
            todo!()
        }
        _ => {
            eprintln!("unknown inlet type '{}'", inlet_type);
            exit(1);
        }
    }
}
