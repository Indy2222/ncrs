use std::convert::Into;

/// Symbol `Other` may represent DNA sequence gaps and misreads.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Symbol {
    Other,
    Adenine,
    Thymine,
    Cytosine,
    Guanine,
}

impl Into<u8> for Symbol {
    fn into(self) -> u8 {
        match self {
            Self::Adenine => 0,
            Self::Thymine => 1,
            Self::Cytosine => 2,
            Self::Guanine => 3,
            Self::Other => 4,
        }
    }
}

impl Into<i32> for Symbol {
    fn into(self) -> i32 {
        let value: u8 = self.into();
        value as i32
    }
}

impl Into<usize> for Symbol {
    fn into(self) -> usize {
        let value: u8 = self.into();
        value as usize
    }
}

/// This struct represents an individual DNA sequencing scaffold, i.e. a
/// continuous sequence of DNA symbols and related metadata.
pub struct Scaffold {
    name: String,
    sequence: Vec<Symbol>,
}

impl Scaffold {
    pub fn new(name: String, sequence: Vec<Symbol>) -> Self {
        Self { name, sequence }
    }

    pub fn name(&self) -> &str {
        self.name.as_str()
    }

    pub fn sequence(&self) -> &[Symbol] {
        &self.sequence
    }
}

/// DNA feature is a human or machine annotated region of a DNA sequence
/// serving a given biological “purpose”. Note that annotations may be mutually
/// overlapping.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Feature {
    Exon,
    /// Protein coding sequence.
    CDS,
    StartCodon,
    StopCodon,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Strand {
    Positive,
    Negative,
}

/// Position of the first symbol (base) of the first full codon/triplet in the
/// feature relative to the feature beginning. Non-zero shift may happen on CDS
/// with start outside of scaffold.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Phase {
    Zero,
    One,
    Two,
}

/// Annotation of a DNA feature.
#[derive(Debug)]
pub struct Annotation {
    scaffold: String,
    source: String,
    feature: Feature,
    score: Option<u32>,
    strand: Strand,
    phase: Option<Phase>,
    start: usize,
    end: usize,
    attributes: String,
}

impl Annotation {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        scaffold: String,
        source: String,
        feature: Feature,
        score: Option<u32>,
        strand: Strand,
        phase: Option<Phase>,
        start: usize,
        end: usize,
        attributes: String,
    ) -> Self {
        Self {
            scaffold,
            source,
            feature,
            score,
            strand,
            phase,
            start,
            end,
            attributes,
        }
    }

    /// Returns identification of the scaffold on which this feature appears.
    pub fn scaffold(&self) -> &str {
        self.scaffold.as_str()
    }

    /// Returns source of the feature (e.g. organization of software name).
    #[allow(dead_code)]
    pub fn source(&self) -> &str {
        self.source.as_str()
    }
    /// Returns type of the annotated feature.
    pub fn feature(&self) -> Feature {
        self.feature
    }

    /// Returns feature quality or confidence.
    #[allow(dead_code)]
    pub fn score(&self) -> Option<u32> {
        self.score
    }

    /// DNA strand on which the feature appears.
    pub fn strand(&self) -> Strand {
        self.strand
    }

    /// See documentation of `Phase`.
    #[allow(dead_code)]
    pub fn phase(&self) -> Option<Phase> {
        self.phase
    }

    /// Inclusive 0-based index of the first symbol of the feature relative to
    /// the scaffold start.
    ///
    /// For example feature `ABC` in sequence `XXABCYYY` would have `.start()`
    /// equal to `2`.
    pub fn start(&self) -> usize {
        self.start
    }

    /// Exclusive 0-based index of the last symbol of the feature relative to
    /// the scaffold start.
    ///
    /// For example feature `ABC` in sequence `XXABCYYY` would have `.end()`
    /// equal to 5.
    pub fn end(&self) -> usize {
        self.end
    }

    /// Attributes of the annotation. Note that the value is take as is and
    /// needs to be further parsed.
    pub fn attributes(&self) -> &str {
        self.attributes.as_str()
    }
}
