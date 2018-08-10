// [[file:~/Workspace/Programming/gchemol/parser/combine.note::b8ed3ff0-660b-4622-89fa-33880a652759][b8ed3ff0-660b-4622-89fa-33880a652759]]
#[macro_use] extern crate combine;
extern crate gchemol;

pub mod xyz;

use combine::{
    char,
    optional,
    Parser,
    Stream
};

use combine::combinator::{
    eof,
    try,
    one_of,
    many,
    many1,
    skip_many,
    skip_many1,
    token,
    take_until,
    recognize
};

use combine::error::{
    ParseError,
};

/// match unsigned integer
pub fn unsigned_integer<I>() -> impl Parser<Input = I, Output = usize>
where
    I: Stream<Item = char>,
    I::Error: ParseError<I::Item, I::Range, I::Position>,
{
    recognize(skip_many1(char::digit())).map(|s: String| {
        s.parse::<usize>().expect("usized integer")
    })
}

/// match signed float number
pub fn signed_float_number<I>() -> impl Parser<Input = I, Output = f64>
where
    I: Stream<Item = char>,
    I::Error: ParseError<I::Item, I::Range, I::Position>,
{
    // +12.6E-4 or - 12.6E-4
    recognize((
        // +12 or -12
        optional(token('+').or(token('-'))),
        skip_many1(char::digit()),
        // .6
        optional((
            token('.'),
            skip_many(char::digit()))),
        // E-4
        optional((
            one_of("Ee".chars()).or(one_of("Dd".chars())),
            optional(token('+').or(token('-'))),
            skip_many1(char::digit()).expected("scientific number")
        ))
    )).map(|s: String| {
        let s = s.replace("D", "E");
        let s = s.replace("d", "e");
        s.parse::<f64>().expect("float number")
    })
}

#[test]
fn test_parser_float() {
    let x = "+124.2D-4";
    let x = signed_float_number().easy_parse(x).unwrap();
    assert_eq!(124.2E-4, x.0);

    let x = "+124.2D+4";
    let x = signed_float_number().easy_parse(x).unwrap();
    assert_eq!(124.2E+4, x.0);

    let x = "-124.2D4";
    let x = signed_float_number().easy_parse(x).unwrap();
    assert_eq!(-124.2E4, x.0);

    let x = "124E-4";
    let x = signed_float_number().easy_parse(x).unwrap();
    assert_eq!(124E-4, x.0);

    let x = "-124.24";
    let x = signed_float_number().easy_parse(x).unwrap();
    assert_eq!(-124.24, x.0);

    let x = "+124.24";
    let x = signed_float_number().easy_parse(x).unwrap();
    assert_eq!(124.24, x.0);

    let x = "+124";
    let x = signed_float_number().easy_parse(x).unwrap();
    assert_eq!(124., x.0);

    let x = "-124.";
    let x = signed_float_number().easy_parse(x).unwrap();
    assert_eq!(-124., x.0);

    let x = "124.";
    let x = signed_float_number().easy_parse(x).unwrap();
    assert_eq!(124., x.0);
}

/// match one or more whitespace char (space or tab)
pub fn ws<I>() -> impl Parser<Input = I, Output = ()>
where
    I: Stream<Item = char>,
    I::Error: ParseError<I::Item, I::Range, I::Position>,
{
    skip_many1(one_of(" \t".chars()))
}

/// match line ending chars: '\n', '\r\n'
pub fn eol<I>() -> impl Parser<Input = I, Output = ()>
where
    I: Stream<Item = char>,
    I::Error: ParseError<I::Item, I::Range, I::Position>,
{
    optional(token('\r')).with(token('\n')).map(|_| ())
    // choice!(
    //     optional(token('\r')).with(token('\n')).map(|_| ()),
    //     eof()
    // )
}

/// read remaing line until line ending (without newline char)
pub fn read_until_eol<I>() -> impl Parser<Input = I, Output = String>
where
    I: Stream<Item = char>,
    I::Error: ParseError<I::Item, I::Range, I::Position>,
{
    (
        take_until(try(eol())),
        eol()
    ).map(|(s, _)| s)
}

#[test]
fn test_parsers() {
    let x = eol().easy_parse("\n").unwrap();
    let x = eol().easy_parse("\r\n").unwrap();
    let x = read_until_eol().easy_parse("wwdd abc\n").unwrap();
}
// b8ed3ff0-660b-4622-89fa-33880a652759 ends here
