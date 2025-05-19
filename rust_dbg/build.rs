// build.rs
use std::process::Command;

fn main() {
    // Get Git commit hash
    let output = Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .unwrap_or_else(|_| panic!("Failed to get Git commit hash"));

    let git_hash = String::from_utf8(output.stdout).unwrap().trim().to_string();

    // Pass it to the compiler
    println!("cargo:rustc-env=GIT_COMMIT_HASH={}", git_hash);
    println!("cargo:rerun-if-changed=.git/HEAD");
}
