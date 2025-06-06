/*! # rust_dbg: a crate for the manipulation of de Bruijn graphs
This crates provides methods for the construction and the manipulation of de Bruijn graphs. It builds heavily upon the crate `debruijn` developped by 10x Genomics.
*/

pub mod fasta_reader;
pub mod graph;

pub mod path;

// pub use ... for re-exports

//####################################################################################
//                              Custom errors                                       //
//####################################################################################
use debruijn::dna_string::DnaString;
use debruijn::Kmer;
use std::error::Error;
use std::io::Write;

/// Custom error type for pathway search operations
#[derive(Debug)]
pub enum PathwayError<K> {
    NoPathExists,
    KmerNotFound(K),
    NodeNotMatching(DnaString, K, usize),
}
impl<K: Kmer> Error for PathwayError<K> {}
impl<K: Kmer> std::fmt::Display for PathwayError<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PathwayError::NoPathExists => write!(f, "No path found between the given k-mers"),
            PathwayError::KmerNotFound(kmer) => write!(f, "Kmer not found: {:?}", kmer),
            PathwayError::NodeNotMatching(seq, kmer, i) => write!(
                f,
                "Expected kmer {:?} at position {} in unitig {}",
                kmer, i, seq
            ),
        }
    }
}

//####################################################################################
//                             Utility functions                                    //
//####################################################################################

/// Parses a string representing a node in the format "(id, direction)"
pub fn parse_node(s: &str) -> Result<(usize, debruijn::Dir), Box<dyn std::error::Error>> {
    let (id, dir) = s
        .strip_prefix('(')
        .and_then(|s| s.strip_suffix(')'))
        .and_then(|s| s.split_once(','))
        .ok_or(format!("Invalid node format: {}", s))?;

    let id = id.parse::<usize>()?;

    let dir = match dir.trim() {
        "Left" => debruijn::Dir::Left,
        "Right" => debruijn::Dir::Right,
        other => return Err(format!("Invalid direction: {}", other).into()),
    };

    Ok((id, dir))
}

/// Removes the previous line and print a new progress bar
pub fn print_progress_bar(current: usize, total: usize) {
    let bar_width = 40;
    let progress = current as f32 / total as f32;
    let filled = (progress * bar_width as f32).round() as usize;
    let empty = bar_width - filled;
    let bar = format!(
        "[{}{}] {}/{} ({:.0}%)",
        "#".repeat(filled),
        " ".repeat(empty),
        current,
        total,
        progress * 100.0
    );
    eprint!("\r{}", bar); // carriage return to overwrite the current line
    std::io::stderr().flush().unwrap();
}

/// Format a long integer with commas
pub fn format_int(n: usize) -> String {
    let s = n.to_string();
    let mut result = String::new();
    let mut chars = s.chars().rev().peekable();

    while let Some(c) = chars.next() {
        result.push(c);
        if chars.peek().is_some() && result.len() % 4 == 3 {
            result.push(',');
        }
    }

    result.chars().rev().collect()
}

//####################################################################################
//                                   Plotting                                       //
//####################################################################################

use plotters::prelude::*;

/// Creates a binned histogram from a vector of usize values and saves it as an SVG file
///
/// # Arguments
/// * `data` - Input data as a vector of usize values
/// * `title` - Title for the histogram plot
/// * `x_min` - Minimum value for x-axis (inclusive)
/// * `x_max` - Maximum value for x-axis (exclusive)
/// * `step` - Bin width
pub fn make_histo(
    data: Vec<usize>,
    title: &str,
    x_min: Option<usize>,
    x_max: Option<usize>,
    step: Option<usize>,
) -> Result<(), Box<dyn Error>> {
    // Determine x_min, x_max, and step
    let x_min = x_min.unwrap_or_else(|| *data.iter().min().unwrap());
    let x_max = x_max.unwrap_or_else(|| *data.iter().max().unwrap() + 1);
    let step = step.unwrap_or(1);

    // Create bins
    let mut bins = Vec::new();
    let mut current = x_min;
    while current < x_max {
        bins.push(current);
        current += step;
    }

    // Initialize counts for each bin
    let mut counts = vec![0; bins.len()];

    // Count values in each bin
    for &value in &data {
        if value >= x_min && value < x_max {
            let bin_index = (value - x_min) / step;
            if bin_index < counts.len() {
                counts[bin_index] += 1;
            }
        }
    }

    // Find maximum count for y-axis
    let y_max = counts.iter().max().copied().unwrap_or(0);

    let file_name = format!("{}.svg", title);
    let drawing_area = SVGBackend::new(&file_name, (600, 400)).into_drawing_area();
    drawing_area.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&drawing_area)
        .caption(title, ("sans-serif", 30))
        .margin(40)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(x_min..x_max, 0..y_max + 1)?;

    chart
        .configure_mesh()
        .x_desc("Value")
        .y_desc("Frequency")
        .axis_desc_style(("sans-serif", 15))
        .draw()?;

    // Draw the histogram bars
    chart.draw_series(bins.iter().zip(counts.iter()).map(|(&x, &y)| {
        Rectangle::new(
            [(x, 0), (x + step, y)], // Each bar spans the bin width
            BLUE.filled(),
        )
    }))?;

    Ok(())
}
