//! reStrainingOrder: mouse strain identification from alignments to an N-masked
//! multi-strain genome. A port of the Perl scripts in the repository root.

pub mod cli;
pub mod order;
pub mod prepare;
pub mod report;
pub mod textio;
pub mod tool;

pub use tool::Tool;
