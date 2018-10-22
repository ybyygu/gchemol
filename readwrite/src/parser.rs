// base

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*base][base:1]]
// line originated parsers
use std::fmt::Debug;

pub use nom::{
    self,
    // Recognizes one or more numerical characters: 0-9
    digit,
    // Recognizes one or more spaces and tabs
    space,
    // Recognizes one or more spaces, tabs, carriage returns and line feeds
    multispace,
    // Recognizes one or more lowercase and uppercase alphabetic characters: a-zA-Z
    alpha,
    alphanumeric,
    is_alphanumeric,
    // Recognizes floating point number in a string and returs a f64
    double as double_s,
    // Recognizes an end of line (both '\n' and '\r\n')
    line_ending,
    not_line_ending,
    // alias
    eol,
};

// Indicating the end of stream
pub const MAGIC_EOF: &str = "$THIS_IS_THE=MAGIC_END_OF_FILE$";
// base:1 ends here

// utils

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*utils][utils:1]]
pub use nom::IResult;
pub use nom::types::{
    CompleteStr,
};

/// whitespace including one or more spaces or tabs
named!(pub space_token<&str, &str>, eat_separator!(&b" \t"[..]));

#[macro_export]
macro_rules! sp (
    ($i:expr, $($args:tt)*) => (
        {
            sep!($i, space_token, $($args)*)
        }
    )
);

/// not any whitespace character
/// will not consume "\n" character
named!(pub not_space<&str, &str>, is_not!(" \t\r\n"));

/// separator using comma or whitespace
named!(pub comma_or_space<&str, &str>, alt!(
    sp!(tag!(",")) | space
));

/// match one unsigned integer
named!(pub unsigned_digit<&str, usize>, map_res!(
    digit,
    str::parse
));

/// match one signed integer
// -1, 0, 1, 2, ...
named!(pub signed_digit<&str, isize>, map_res!(
    recognize!(
        pair!(
            opt!(alt!(char!('+') | char!('-'))),
            digit
        )
    ),
    str::parse
));

// sections are separated by a blank line
named!(pub blank_line<&str, &str>, sp!(line_ending));

#[test]
fn test_parser_signed_digit() {
    let (_, x) = signed_digit("12\n")
        .expect("parser: signed_digit 12");
    assert_eq!(12, x);

    let (_, x) = signed_digit("+12\n")
        .expect("parser: signed_digit +12");
    assert_eq!(12, x);

    let (_, x) = signed_digit("-12\n")
        .expect("parser: signed_digit -12");
    assert_eq!(-12, x);
}

/// match one or more unsigned numbers separated by whitespace
named!(pub count_many<&str, Vec<usize>>, terminated!(
    many1!(sp!(unsigned_digit)),
    sp!(eol)
));

#[test]
fn test_parser_count_many() {
    let (_, ns) = count_many(" 1 2 3 4 5 \n")
        .expect("parser: count_many");
    assert_eq!(5, ns.len());
}

/// read the remaining line including the eol character
named!(pub read_until_eol<&str, &str>, terminated!(
    not_line_ending,
    eol
));

#[test]
fn test_parser_read_until_eol() {
    let x = read_until_eol("this is the end\nok\n")
        .expect("parser: read_until_eol");
    let x = read_until_eol("\n")
        .expect("parser: read_until_eol empty line");
}

/// Consume three float numbers separated by one or more spaces
/// Return position array
named!(pub xyz_array<&str, [f64; 3]>, do_parse!(
    x: double_s >>
       space    >>
    y: double_s >>
       space    >>
    z: double_s >>
    (
        [x, y, z]
    )
));

#[test]
fn test_parser_xyz_array() {
    let (_, x) = xyz_array("-11.4286  1.7645  0.0000 ").unwrap();
    assert_eq!(x[2], 0.0);

    let (_, x) = xyz_array("-11.4286  1.7645  0.0000\n").unwrap();
    assert_eq!(x[2], 0.0);

    let (_, x) = xyz_array("-11.4286\t1.7E-5  0.0000 \n").unwrap();
    assert_eq!(x[2], 0.0);
}
// utils:1 ends here

// eof

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*eof][eof:1]]
named!(pub magic_eof<&str, &str>, tag!(MAGIC_EOF));

// #[test]
// fn test_magic_eof() {
//     let x = "aaon";
//     let x = eof(x);
//     println!("{:#?}", x);
// }
// eof:1 ends here
