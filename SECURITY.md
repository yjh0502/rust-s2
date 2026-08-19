# Security Policy 🔒

Thanks for helping keep `s2` and its users safe! If you've found a security
vulnerability, we'd genuinely like to hear about it.

## Reporting a vulnerability 📬

Please report security issues **privately**, by emailing
**sascha@brawer.ch**, rather than opening a public issue. Include as much
detail as you can — affected version(s), a reproduction if possible, and
the potential impact — so it's easier to reproduce and fix quickly.

(GitHub also offers a built-in [private vulnerability
reporting](https://docs.github.com/en/code-security/security-advisories/guidance-on-reporting-and-writing-information-about-vulnerabilities/privately-reporting-a-security-vulnerability)
flow via a repository's Security tab. It isn't usable here yet — reports
made through it aren't currently visible to the active maintainer — so
please use email instead until that changes.)

This is a small, mostly-one-person-maintained project 🙂 — there's no
formal SLA, but security reports get priority attention over everything
else in the queue.

## Supported versions

Only the most recently published version on
[crates.io](https://crates.io/crates/s2) is supported. If you're on an
older version, please upgrade before reporting — the issue may already be
fixed.

## Supply-chain & release security 🔗

A few things worth knowing about how `s2` gets built and published, for
anyone doing their own due diligence:

- **Dependency scanning**: every dependency is checked on every PR and
  weekly (`.github/workflows/cargo-deny.yml`) for known security
  advisories, incompatible licenses, and untrusted sources — see
  `deny.toml`.
- **Trusted Publishing**: releases are published to crates.io via
  [OIDC-based Trusted
  Publishing](https://crates.io/docs/trusted-publishing) rather than a
  long-lived API token — there's no static publishing credential that
  could leak or be stolen. A release can only be published by the exact,
  pinned GitHub Actions workflow in this repository.
- **SLSA Build Level 3 provenance**: every release package is
  [attested](https://slsa.dev/spec/v1.2/build-track-basics#build-l3) —
  signed, tamper-evident proof of exactly which commit and workflow
  produced it, generated in a build environment isolated from the actual
  build/package steps so a compromised build can't forge its own
  attestation. Verify a downloaded `.crate` file yourself with
  [`gh attestation
  verify`](https://cli.github.com/manual/gh_attestation_verify)
  (`gh attestation verify s2-X.Y.Z.crate --repo yjh0502/rust-s2`). Note
  that crates.io itself doesn't check this automatically — `cargo
  add`/`cargo build` won't verify it for you; this is for anyone who
  wants to check for themselves.

For how to report a *non-security* bug, see
[CONTRIBUTING.md](CONTRIBUTING.md) instead.
