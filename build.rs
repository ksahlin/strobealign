fn main() {
    cc::Build::new()
        .file("ext/ssw/ssw.c")
        .compile("ssw");
}
