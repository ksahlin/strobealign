use std::error::Error;
use std::fs::read;
use std::io::Write;

use assert_cmd::{self, Command};
use flate2::Compression;
use flate2::write::GzEncoder;
use predicates::prelude::*;
use temp_file::TempFileBuilder;

fn cmd() -> Command {
    Command::new(assert_cmd::cargo::cargo_bin!())
}

#[test]
fn fail_without_arguments() {
    let mut cmd = cmd();
    cmd.assert().failure();
}

#[test]
fn fail_with_unknown_argument() {
    let mut cmd = cmd();
    cmd.arg("-G").assert().failure();
}

#[test]
fn success_when_printing_help() {
    let mut cmd = cmd();
    cmd.arg("--help").assert().success();
    cmd.arg("-h").assert().success();
}

/// Without --eqx, there must be no = and X CIGAR operators
#[test]
fn single_end_sam_m_cigar() {
    let mut cmd = cmd();
    let p = cmd.args(&["tests/phix.fasta", "tests/phix.1.fastq"]);
    let a = p.assert().success();
    let stdout = &a.get_output().stdout;
    for line in str::from_utf8(stdout).unwrap().lines() {
        if line.as_bytes()[0] == b'@' {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let cigar = fields[5];
        assert!(cigar.find('=').is_none());
        assert!(cigar.find('X').is_none());
    }
}

#[test]
fn single_end_sam() {
    // TODO -v
    let expected = String::from_utf8(read("tests/phix.se.sam").unwrap()).unwrap();
    let mut cmd = cmd();
    let p = cmd.args(&[
        "--chunk-size=3",
        "--no-PG",
        "--eqx",
        "--rg-id=1",
        "--rg=SM:sample",
        "--rg=LB:library",
        "tests/phix.fasta",
        "tests/phix.1.fastq",
    ]);
    let a = p.assert();
    a.success().stdout(predicate::str::diff(expected));
}

#[test]
fn paired_end_sam() {
    let expected = String::from_utf8(read("tests/phix.pe.sam").unwrap()).unwrap();
    let mut cmd = cmd();
    let p = cmd.args(&[
        "--no-PG",
        "--eqx",
        "--rg-id=1",
        "--rg=SM:sample",
        "--rg=LB:library",
        "tests/phix.fasta",
        "tests/phix.1.fastq",
        "tests/phix.2.fastq",
    ]);
    let a = p.assert();
    a.success().stdout(predicate::str::diff(expected));
}

#[test]
fn compressed_reference() -> Result<(), Box<dyn Error>> {
    let mut cmd = cmd();
    let fasta_contents = read("tests/phix.fasta")?;
    let r = String::from_utf8(fasta_contents.clone())?;
    dbg!(r);
    let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
    encoder.write_all(&fasta_contents)?;
    let compressed_fasta = encoder.finish()?;
    dbg!(&compressed_fasta);
    let tmp = TempFileBuilder::new()
        .suffix(".gz")
        .build()?
        .with_contents(&compressed_fasta)?;

    cmd.args(&[
        tmp.path().as_os_str().to_str().unwrap(),
        "tests/empty.fastq",
    ])
    .assert()
    .success();

    Ok(())
}
