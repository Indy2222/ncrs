use crate::data::{Annotation, Feature, Phase, Strand};
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::path::Path;

const GFF_NUM_COLUMNS: usize = 9;

/// Load scaffold annotations from a general feature format (GFF) file.
pub fn load_gff_file(path: &Path) -> Result<Vec<Annotation>> {
    let reader = {
        let file =
            File::open(path).with_context(|| format!("Could not open file {}.", path.display()))?;
        BufReader::new(file)
    };

    let mut annotations = Vec::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("Could not read file {}.", path.display()))?;
        let annotation = parse_gff_line(line).with_context(|| {
            format!("Failed to parse line {} of file {}.", i + 1, path.display())
        })?;
        annotations.push(annotation);
    }

    Ok(annotations)
}

fn parse_gff_line(line: String) -> Result<Annotation> {
    let mut tokens: Vec<String> = line.split('\t').take(9).map(String::from).collect();

    let num_columns = tokens.len();
    if num_columns != GFF_NUM_COLUMNS {
        return Err(anyhow!(
            "Not enough tab separated tokens. Expected {} got {}.",
            GFF_NUM_COLUMNS,
            num_columns
        ));
    }

    let attributes = tokens.pop().unwrap();
    let phase = tokens.pop().unwrap();
    let strand = tokens.pop().unwrap();
    let score = tokens.pop().unwrap();
    let end = tokens.pop().unwrap();
    let start = tokens.pop().unwrap();
    let feature = tokens.pop().unwrap();
    let source = tokens.pop().unwrap();
    let scaffold = tokens.pop().unwrap();
    assert_eq!(tokens.len(), 0);

    let phase = match phase.as_str() {
        "0" => Some(Phase::Zero),
        "1" => Some(Phase::One),
        "2" => Some(Phase::Two),
        _ => None,
    };

    let strand = match strand.as_str() {
        "+" => Strand::Positive,
        "-" => Strand::Negative,
        unrecognized => {
            return Err(anyhow!(
                "Invalid strand, only +, - are valid. Got: {}",
                unrecognized
            ));
        }
    };

    let score = match score.as_str() {
        "." => None,
        score => Some(
            score
                .parse::<u32>()
                .context("Score is not a positive integer.")?,
        ),
    };

    // GFF end is 1-based inclusive, we want 0-based exclusive which is the same number.
    let end = end
        .parse()
        .with_context(|| format!("Feature end has to be a positive integer. Got: {}", end))?;

    let start = start.parse::<usize>().with_context(|| {
        format!(
            "Annotation start has to be a positive integer. Got: {}",
            start
        )
    })?;

    if start == 0 {
        return Err(anyhow!(
            "Feature start position is 0 but must be bigger or equal to 1."
        ));
    }
    // GFF start is 1-based inclusive and we want 0-based inclusive, therefore
    // the subtraction.
    let start = start - 1;

    if start >= end {
        return Err(anyhow!(
            "Feature start index is greater or equal to end index. {} >= {}",
            start,
            end
        ));
    }

    let feature = match feature.as_str() {
        "start_codon" => Feature::StartCodon,
        "stop_codon" => Feature::StopCodon,
        "CDS" => Feature::CDS,
        "exon" => Feature::Exon,
        unrecognized => {
            return Err(anyhow!("Unrecognized feature: {}", unrecognized));
        }
    };

    Ok(Annotation::new(
        scaffold, source, feature, score, strand, phase, start, end, attributes,
    ))
}

#[cfg(test)]
mod test {

    use crate::data::{Feature, Phase, Strand};
    use std::path::Path;

    #[test]
    fn test_load_valid_gff() {
        let gff_path = Path::new("./tests/valid.gff");
        let mut annotations = super::load_gff_file(gff_path).unwrap();
        assert_eq!(annotations.len(), 4);

        let four = annotations.pop().unwrap();
        let three = annotations.pop().unwrap();
        let two = annotations.pop().unwrap();
        let one = annotations.pop().unwrap();

        assert_eq!(one.scaffold(), "scaffold_1");
        assert_eq!(one.feature(), Feature::Exon);

        assert_eq!(two.scaffold(), "scaffold_2");
        assert_eq!(two.feature(), Feature::CDS);

        assert_eq!(three.scaffold(), "scaffold_3");
        assert_eq!(three.feature(), Feature::StartCodon);

        assert_eq!(four.scaffold(), "scaffold_4");
        assert_eq!(four.source(), "JGI");
        assert_eq!(four.feature(), Feature::StopCodon);
        assert_eq!(four.start(), 2183);
        assert_eq!(four.end(), 2186);
        assert_eq!(four.score(), None);
        assert_eq!(four.strand(), Strand::Positive);
        assert_eq!(four.phase(), Some(Phase::Zero));

        assert_eq!(
            four.attributes(),
            "name \"fgenesh1_kg.1_#_1_#_Locus4417v1rpkm26.65\""
        );
    }

    #[test]
    fn test_load_invalid_gff() {
        let gff_path = Path::new("./tests/invalid.gff");
        let error = match super::load_gff_file(gff_path) {
            Ok(_) => panic!("Loading did not fail."),
            Err(error) => error,
        };

        assert_eq!(
            format!("{}", error),
            String::from("Failed to parse line 2 of file ./tests/invalid.gff.")
        );

        assert_eq!(
            format!("{}", error.root_cause()),
            String::from("Unrecognized feature: XXX")
        );
    }
}
