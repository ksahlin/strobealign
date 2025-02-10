use assert_cmd::Command;
use std::fs::read;
use predicates::prelude::*;

#[test]
fn fail_without_arguments() {
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.assert().failure();
}

#[test]
fn fail_with_unknown_argument() {
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.arg("-G").assert().failure();
}

#[test]
fn success_when_printing_help() {
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.arg("--help").assert().success();
    cmd.arg("-h").assert().success();
}

#[test]
fn single_end_sam() {
    // TODO --chunk-size 3
    // TODO -v
    let expected = String::from_utf8(read("tests/phix.se.sam").unwrap()).unwrap();
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    let p = cmd.args(&["--no-PG", "--eqx", "--rg-id=1", "--rg=SM:sample", "--rg=LB:library", "tests/phix.fasta", "tests/phix.1.fastq"]);
    let a = p.assert();
    a.success().stdout(predicate::str::diff(expected));
}

#[test]
fn paired_end_sam() {
    // TODO --chunk-size 3
    let expected = String::from_utf8(read("tests/phix.pe.sam").unwrap()).unwrap();
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    let p = cmd.args(&["--no-PG", "--eqx", "--rg-id=1", "--rg=SM:sample", "--rg=LB:library", "tests/phix.fasta", "tests/phix.1.fastq", "tests/phix.2.fastq"]);
    let a = p.assert();
    a.success().stdout(predicate::str::diff(expected));
}