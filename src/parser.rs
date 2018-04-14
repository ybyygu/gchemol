// [[file:~/Workspace/Programming/gchemol/gchemol.note::85054519-f2d5-4c63-994a-78bbe4f9a30f][85054519-f2d5-4c63-994a-78bbe4f9a30f]]
use std::fmt::Debug;
use nom::IResult;
use nom::{alphanumeric, multispace, float_s, double_s, is_digit, digit, crlf, line_ending};

named!(maybe_f64_s<&str, f64>,
       map_res!(is_not_s!(" \t"), str::parse)
);

named!(pub xyz_array<&str, [f64; 3]>,
       do_parse!(
           x: ws!(maybe_f64_s) >>
           y: ws!(maybe_f64_s) >>
           z: ws!(maybe_f64_s) >>
           ([x, y, z])
       )
);

named!(usize_timestep_s<&str, usize>,
       ws!(preceded!(
           tag!("# Timestep"),
           map_res!(digit, str::parse)
       ))
);

named!(usize_nparticles_s<&str, usize>,
       ws!(preceded!(
           tag!("# Number of particles"),
           map_res!(digit, str::parse)
       ))
);

named!(usize_neighbors_s<&str, Vec<usize>>,
       length_count!(
           usize_timestep_s,
           map_res!(ws!(digit), str::parse)
       )
);

/// return last usize digit in a line
/// 200\n => 200_usize
named!(usize_last_digit<&str, usize>,
       map_res!(
           terminated!(
               digit,
               alt!(eof!() | line_ending)
           ),
           str::parse
       )
);

named!(last_item<&str, &str>,
       terminated!(
           is_not_s!(" \t\n\r"),
           alt!(line_ending | eof!())
       )
);

// # Timestep 3 => 3
named!(pub take_last_digit<&str, usize>,
       do_parse!(
           r: many_till!(take!(1), ws!(usize_last_digit)) >>
           (r.1)
       )
);

// # Timestep 3 => 3
named!(pub take_last_item<&str, &str>,
       do_parse!(
           r: many_till!(take!(1), ws!(last_item)) >>
           (r.1)
       )
);
// 85054519-f2d5-4c63-994a-78bbe4f9a30f ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::673f851b-fd18-4dc1-9598-80e6dc83cf14][673f851b-fd18-4dc1-9598-80e6dc83cf14]]
#[test]
fn test_nom() {
    let x = usize_timestep_s("# Timestep 50");
    println!("{:?}", x);
    let x = usize_neighbors_s(" 408 1 8 416 272 280 288 400 407 392 536 0 0.803 0.059");
    let x = take_last_digit("# Number of particles 50");
    let x = take_last_item("# Number of particles 50.2
  ");

    let (_, x) = xyz_array(" 1. 1. 2").unwrap();
    assert_eq!([1.0, 1.0, 2.0], x);
}
// 673f851b-fd18-4dc1-9598-80e6dc83cf14 ends here
