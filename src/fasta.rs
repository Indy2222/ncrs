use crate::data::{Scaffold, Symbol};
use anyhow::{Context, Result};

use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::path::Path;

struct ScaffoldBuilder {
    name: String,
    sequence: Vec<Symbol>,
}

impl ScaffoldBuilder {
    fn new(name: String) -> Self {
        Self {
            name,
            sequence: Vec::new(),
        }
    }

    fn extend_from_str(&mut self, seq: &str) -> Result<()> {
        let seq = seq
            .chars()
            .map(|c| match c {
                'A' | 'a' => Ok(Symbol::Adenine),
                'C' | 'c' => Ok(Symbol::Cytosine),
                'T' | 't' => Ok(Symbol::Thymine),
                'G' | 'g' => Ok(Symbol::Guanine),
                'N' | 'n' => Ok(Symbol::Other),
                _ => Err(anyhow!("Encountered invalid symbol {}.", c)),
            })
            .collect::<Result<Vec<Symbol>>>()?;
        self.sequence.extend(seq);
        Ok(())
    }

    fn build(self) -> Scaffold {
        let Self { name, sequence } = self;
        Scaffold::new(name, sequence)
    }
}

/// Load FASTA file.
pub fn load_fasta(path: &Path) -> Result<Vec<Scaffold>> {
    let mut reader = {
        let file =
            File::open(path).with_context(|| format!("Failed to open file {}.", path.display()))?;
        BufReader::new(file)
    };

    let mut scaffolds = Vec::new();

    let mut builder: Option<ScaffoldBuilder> = None;
    let mut line = String::new();

    loop {
        let num_bytes = reader
            .read_line(&mut line)
            .with_context(|| format!("Failed to read file {}.", path.display()))?;

        if num_bytes == 0 {
            break;
        }

        if line.ends_with('\n') {
            line.pop();
            if line.ends_with('\r') {
                line.pop();
            }
        }

        if line.starts_with('>') {
            if let Some(builder) = builder {
                scaffolds.push(builder.build());
            }

            builder = Some(ScaffoldBuilder::new(String::from(&line[1..])));
        } else {
            match builder {
                Some(ref mut b) => b.extend_from_str(&line)?,
                None => {
                    return Err(anyhow!("Ivalid FASTA file {}.", path.display()));
                }
            }
        }

        line.clear();
    }

    match builder {
        Some(builder) => scaffolds.push(builder.build()),
        None => return Err(anyhow!("Empty FASTA file {}.", path.display())),
    }

    Ok(scaffolds)
}

#[cfg(test)]
mod test {

    use crate::data::Symbol;
    use std::path::Path;

    #[test]
    fn test_load_fasta() {
        let fasta_path = Path::new("./tests/valid.fasta");
        let mut scaffolds = super::load_fasta(fasta_path).unwrap();
        assert_eq!(scaffolds.len(), 2);

        let second = scaffolds.pop().unwrap();
        let first = scaffolds.pop().unwrap();

        assert_eq!(first.name(), "scaffold_1");
        assert_eq!(first.sequence().len(), 280);
        assert_eq!(second.name(), "scaffold_2");

        let expected_sequence = vec![
            Symbol::Thymine,
            Symbol::Thymine,
            Symbol::Cytosine,
            Symbol::Thymine,
            Symbol::Guanine,
            Symbol::Other,
            Symbol::Adenine,
        ];
        assert_eq!(second.sequence(), &expected_sequence[..]);
    }
}
