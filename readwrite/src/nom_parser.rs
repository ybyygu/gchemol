// base

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*base][base:1]]
// Indicating the end of stream
pub const MAGIC_EOF: &str = "\n\nxTHIS_IS_THE=MAGIC_END_OF_FILE\n";
// base:1 ends here

// separator

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*separator][separator:1]]
named!(pub eof<&str, &str>, tag!(MAGIC_EOF));

/// A whitespace wrapper consuming " \t\r" (no newline)
named!(pub space_token<&str, &str>, eat_separator!(&b" \t\r"[..]));

#[macro_export]
macro_rules! sp (
    ($i:expr, $($args:tt)*) => (
        {
            sep!($i, space_token, $($args)*)
        }
    )
);

/// A whitespace wrapper consuming "\r\n" (line-ending)
named!(pub eol_token<&str, &str>, eat_separator!(&b"\r\n"[..]));

#[macro_export]
macro_rules! nl (
    ($i:expr, $($args:tt)*) => (
        {
            sep!($i, eol_token, $($args)*)
        }
    )
);

/// A DOS `\r` char consumer
named!(pub dos_eol_pre_token<&str, &str>, eat_separator!(&b"\r\n"[..]));

#[macro_export]
macro_rules! dos2unix (
    ($i:expr, $($args:tt)*) => (
        {
            sep!($i, dos_eol_pre_token, $($args)*)
        }
    )
);

#[test]
fn test_parser_sp() {
    // will consume space between each token.
    let _  = sp!(" x a b c\t x d",
                tuple!(tag!("x"),
                       tag!("a"),
                       tag!("b"),
                       tag!("c")
                )).expect("parser sp test1");


    // automatically remove line carriage char `\r` in DOS format
    let _ = sp!(" x\t a b c\r\nx d\r\n",
                     pair!(
                         tuple!(tag!("x"),
                                tag!("a"),
                                tag!("b"),
                                tag!("c")
                         ),
                         tag!("\nx")
                     )
    ).expect("parser sp test2");
}
// separator:1 ends here

// error
// How to use:
// : return Err(nom_error!(input));
// : return Err(nom_error!(input, 120));
// : map_err(|e| nom_failure!(input))

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*error][error:1]]
#[macro_export]
macro_rules! nom_error {
    ($e:expr) => {
        nom::Err::Error(error_position!($e, nom::ErrorKind::Custom(28)))
    };
    ($e:expr, $c:expr) => {
        nom::Err::Error(error_position!($e, nom::ErrorKind::Custom($c)))
    };
}

#[macro_export]
macro_rules! nom_failure {
    ($e:expr) => {
        nom::Err::Failure(error_position!($e, nom::ErrorKind::Custom(38)))
    };
    ($e:expr, $c:expr) => {
        nom::Err::Failure(error_position!($e, nom::ErrorKind::Custom($c)))
    };
}
// error:1 ends here

// reexport
// reexport nom combinators


// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*reexport][reexport:1]]
pub use nom::{
    self,
    // Recognizes floating point number in a string and returs a f64
    double,
    // Recognizes one or more numerical characters: 0-9
    digit,
    // Recognizes one or more spaces and tabs
    space,
    // Recognizes one or more spaces, tabs, carriage returns and line feeds
    multispace,
    // Recognizes one or more lowercase and uppercase alphabetic characters: a-zA-Z
    alpha,
    alphanumeric,
    // Recognizes an end of line (both '\n' and '\r\n')
    line_ending,
    // Shorter alias
    eol,
    // Everything except eol
    not_line_ending,
};
// reexport:1 ends here

// separators

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*separators][separators:1]]
// Match a blank line containing zero or more whitespace character
named!(pub blank_line<&str, &str>, sp!(nom::line_ending));

/// Anything except whitespace
/// will not consume "\n" character
named!(pub not_space<&str, &str>, is_not!(" \t\r\n"));

/// separator using comma or whitespace
named!(pub comma_or_space<&str, &str>, alt!(
    sp!(tag!(",")) | space
));
// separators:1 ends here

// numbers

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*numbers][numbers:1]]
/// Match one unsigned integer: 123
named!(pub unsigned_digit<&str, usize>, map_res!(
    digit,
    str::parse
));

/// match one signed integer: -1, 0, 1, 2, ...
named!(pub signed_digit<&str, isize>, map_res!(
    recognize!(
        pair!(
            opt!(alt!(char!('+') | char!('-'))),
            digit
        )
    ),
    str::parse
));

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
// numbers:1 ends here

// coordinates

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*coordinates][coordinates:1]]
/// Consume three float numbers separated by one or more spaces
/// Return position array
named!(pub xyz_array<&str, [f64; 3]>, do_parse!(
    x: double   >>
       space    >>
    y: double   >>
       space    >>
    z: double   >>
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
// coordinates:1 ends here

// numbers

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*numbers][numbers:1]]
/// Parse a line containing an unsigned integer
named!(pub read_usize<&str, usize>, sp!(terminated!(unsigned_digit, eol)));

/// Parse a line containing many unsigned integers
named!(pub read_usize_many<&str, Vec<usize>>,
       sp!(terminated!(
           many1!(unsigned_digit),
           eol
       ))
);

#[test]
fn test_parser_usize_many() {
    let (_, ns) = read_usize_many(" 11 2 3 4 5 \r\n")
        .expect("parser: count_many");
    assert_eq!(5, ns.len());
}

/// Parse a line containing a float number
named!(pub read_f64<&str, f64>, sp!(terminated!(double, eol)));

/// Parse a line containing many float numbers
named!(pub read_f64_many<&str, Vec<f64>>, sp!(
    terminated!(
        many1!(double),
        eol
    )
));

#[test]
fn test_parser_f64_many() {
    let line = "1.2  3.4 -5.7 0.2 \n";
    let (_, fs) = read_f64_many(line).expect("f64 parser");
    assert_eq!(4, fs.len());
}
// numbers:1 ends here

// lines

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*lines][lines:1]]
/// Match the remaining line including the eol (end of line) character
named!(pub read_until_eol<&str, &str>, terminated!(
    not_line_ending,
    line_ending
));

/// Peek current line without consuming it
named!(pub peek_line<&str, &str>, peek!(terminated!(
    not_line_ending,
    line_ending
)));

// /// # Note
// /// Use `sp!` macro to remove the remaining `\r` char in DOS format
// named!(pub read_line<&str, &str>, take_until_and_consume!("\n"));

/// Read a single line, possible faster than read_until_eol?
#[inline]
pub fn read_line(input: &str) -> nom::IResult<&str, &str> {
    let (rest, line) = take_until_and_consume!(input, "\n")?;

    // remove remaining carriage return: `\r` for windows line ending
    let n = line.len();
    if line.ends_with("\r") {
        Ok((rest, &line[0..n-1]))
    } else {
        Ok((rest, line))
    }
}

#[test]
fn test_parser_read_until_eol() {
    let _ = read_until_eol("this is the end\nok\n").expect("parser: read_until_eol");
    let _ = read_until_eol("\n").expect("parser: read_until_eol empty line");

    let (rest, line) = read_line("this is the end\r\nok\r\n")
        .expect("parser: read_until_eol");
    assert_eq!("this is the end", line);
    assert_eq!("ok\r\n", rest);

    let (rest, line) = read_line("\n\n")
        .expect("parser: read_line empty line");
    assert_eq!("", line);
    assert_eq!("\n", rest);
}

/// Read m lines from the stream
#[inline]
pub fn read_many_lines(input: &str, m: usize) -> nom::IResult<&str, Vec<&str>> {
    many_m_n!(input, m, m, read_line)
}

#[test]
fn test_parser_read_many_lines() {
    let txt = "12
test
C -11.4286  1.7645  0.0000
C -10.0949  0.9945  0.0000
C -10.0949 -0.5455  0.0000
C -11.4286 -1.3155  0.0000
";
    let (_, lines) = read_many_lines(txt, 3).expect("read_many_lines");
    assert_eq!(3, lines.len());
}

/// Read lines until the line starting with the `label`, and return the consumed lines
#[inline]
pub fn read_lines_until<'a>(input: &'a str, label: &'a str) -> nom::IResult<&'a str, Vec<&'a str>> {
    let (input, pp) = dos2unix!(
        input,
        many_till!(
            read_line,
            peek!(
                alt!(tag!(label) | eof)
            )
        )
    )?;

    Ok((input, pp.0))
}

/// Read lines until the line starting with the `label` and consumes it
#[inline]
pub fn read_lines_until_and_consume<'a>(input: &'a str, label: &'a str) -> nom::IResult<&'a str, Vec<&'a str>> {
    let (input, pp) = dos2unix!(input,
        many_till!(
            read_line,
            alt!(tag!(label) | eof)
        )
    )?;

    Ok((input, pp.0))
}

/// Read lines until the parser applies
#[inline]
pub fn read_many_until<'a, F>(input: &'a str, parser: F) -> nom::IResult<&'a str, Vec<&'a str>>
where
    F: Fn(&'a str) -> nom::IResult<&'a str, &'a str>,
{
    let (input, pp) = many_till!(
        input,
        read_line,
        peek!(
            alt!(
                parser |
                eof
            )
        )
    )?;

    Ok((input, pp.0))
}

#[inline]
pub fn read_many_until_and_consume<'a, F>(input: &'a str, parser: F) -> nom::IResult<&'a str, Vec<&'a str>>
where
    F: Fn(&'a str) -> nom::IResult<&'a str, &'a str>,
{
    let (input, pp) = many_till!(
        input,
        read_line,
        alt!(
            parser |
            eof
        )
    )?;

    Ok((input, pp.0))
}
// lines:1 ends here
