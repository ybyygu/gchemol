// [[file:~/Workspace/Programming/gchemol/parser/combine.note::b8ed3ff0-660b-4622-89fa-33880a652759][b8ed3ff0-660b-4622-89fa-33880a652759]]
#[macro_use] extern crate combine;

use combine::{
    char,
    optional,
    choice,
    Parser,
    Stream
};

use combine::combinator::{
    one_of,
    value,
    satisfy,
    skip_many,
    skip_many1,
    token,
    recognize
};

use combine::error::{
    ParseError,
    UnexpectedParse,
};

fn signed_float_number<I>() -> impl Parser<Input = I, Output = f64>
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

    let x = "124E-4";
    let x = signed_float_number().easy_parse(x).unwrap();
    assert_eq!(124E-4, x.0);
}
// b8ed3ff0-660b-4622-89fa-33880a652759 ends here
