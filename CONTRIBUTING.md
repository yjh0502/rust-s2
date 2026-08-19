# Contributing to rust-s2 🌐

Thanks for taking a look — contributions of any size are genuinely welcome
here, whether that’s a typo fix, a bug report, a new test, or porting an
entire piece of functionality. This document exists so you don’t have to
guess at the conventions; skim the parts that are relevant to what you’re
doing and don’t worry about the rest.

## Background & guiding principles 🧭

`s2` traces back through two ancestors:

- [**google/s2geometry**](https://github.com/google/s2geometry) — the
  original C++ library, and the authoritative source of truth for what the
  S2 algorithms actually do and why.
- [**golang/geo**](https://github.com/golang/geo) — a Go port of the C++
  library, which is what this crate was originally ported *from*.

Looking at the Go port is a fine place to start, but when you’re porting
or fixing something, please also check the **C++ source** — it’s the
authoritative implementation, generally better maintained than the Go
port, which is a step removed from it and occasionally trails or diverges
in subtle ways. If you’re porting a function, port its tests too (adapted
to Rust, not copy-pasted): the C++ and Go test suites encode a lot of
hard-won edge-case knowledge (degenerate geometry, numerical precision
boundaries, antipodal points, and the like) that’s easy to miss if you’re
only working from the algorithm description.

That said, this is **not meant to be a literal, line-by-line port**. The
goal is a Rust API that feels natural to a Rust programmer — using
`Result`/`Option`, iterators, trait implementations, and ownership the way
idiomatic Rust code would, rather than mirroring the C++/Go class and
method structure just because that’s how the source happened to be
organized. When the idiomatic Rust shape and the original structure
diverge, prefer the idiomatic Rust shape, and leave a short comment
pointing at what you ported from so the next person can find their way
back to the source of truth.

## Getting started 🚀

1. Fork the repo and create a branch for your change.
2. Make your change (see the guiding principles above if you’re porting or
   fixing S2 algorithm code).
3. Before opening a PR, run what CI will check anyway:
   ```sh
   cargo fmt --all
   cargo clippy --all-features --all-targets
   cargo test --all-features
   ```
4. Open a pull request — see the title convention below. Everything else
   CI will guide you through.

New to the codebase? Issues labeled [`good first
issue`](https://github.com/yjh0502/rust-s2/labels/good%20first%20issue)
are a good place to start. Feel free to open an issue first if you’d like
to talk through an approach before writing code — nobody will mind.

CI also runs two checks that aren’t in the list above, so don’t be
surprised if either trips you up even though everything passed locally:

- **`cargo-deny`** (`deny.toml`): new dependencies must come from
  crates.io, use a permissively-licensed license already on the allow
  list (currently MIT/Apache-2.0/Unicode-3.0), and have no known
  security advisory. Install with `cargo install cargo-deny` and run
  `cargo deny check` to verify ahead of time.
- **MSRV build**: CI also compiles on the exact minimum Rust version this
  crate declares (`rust-version` in `Cargo.toml`, currently 1.86) — a
  standard-library API newer than that will build fine on your machine
  and still fail this job.

## Pull request titles: Conventional Commits 📝

This repo squash-merges PRs, so your PR title becomes the commit message
on `master`. Please title PRs using [Conventional
Commits](https://www.conventionalcommits.org/en/v1.0.0/):

```
<type>[optional scope][!]: <description>
```

- `feat:` — a new, backward-compatible public API addition or capability.
- `fix:` — a bug fix.
- `chore:`, `ci:`, `refactor:`, `test:`, `build:`, `docs:` — everything with
  no public-facing effect: CI/tooling changes, internal refactors, tests,
  documentation. These don’t trigger a release version bump.
- Append `!` right after the type (and optional scope) for anything that
  breaks the public API or raises the
  [MSRV](https://doc.rust-lang.org/cargo/reference/rust-version.html), e.g.
  `feat!:`, `fix!:`, `chore!:`. Breaking changes bump the *minor* version
  while this crate is pre-1.0 (`0.x.y` → `0.(x+1).0`), per [Cargo’s SemVer
  compatibility](https://doc.rust-lang.org/cargo/reference/semver.html) —
  not the (nonexistent, for a `0.x` crate) major version.

Don’t worry about getting this perfect on the first try: CI checks it
(`.github/workflows/pr-title-lint.yml`) and will tell you if the title
needs adjusting, right on the PR.

This convention exists to eventually feed
[release-plz](https://release-plz.dev), which can fully automate version
bumps, `CHANGELOG.md` generation, tagging, and publishing from commits
formatted this way. It isn’t wired up yet, pending a possible repo
ownership move (see [#49](https://github.com/yjh0502/rust-s2/issues/49)).
Releases are cut manually until then (see below), but PR titles should
already follow this convention so that enabling release-plz later is a
small activation step, not a retrofit.

## GitHub labels 🏷️

PRs are also labeled for GitHub’s auto-generated release notes
(`.github/release.yml`). CI derives and applies the label from your PR
title automatically (`.github/workflows/pr-title-lint.yml`), so in the
common case there’s nothing to do here either:

| Label | Meaning | Auto-applied from |
|---|---|---|
| `breaking-change` | Public API changed incompatibly | `!` on the commit type |
| `enhancement` | New public API/capability, non-breaking | `feat:` |
| `bug` | Bug fix | `fix:` |
| `msrv` | Minimum supported Rust version raised | *(manual — see below)* |
| *(none)* | Purely internal, no public-facing effect | `chore:`/`ci:`/`refactor:`/`test:`/`build:`/`docs:` |

The one exception is `msrv`: it’s a narrow, infrequent special case (an
MSRV bump is also always `chore!:` and thus already gets `breaking-change`
automatically) that only a human can reliably recognize, so please add it
by hand on top when it applies.

### How to tell whether a change is actually breaking

A diff can look identical whether or not the item it touches was ever
reachable from outside the crate — so don’t infer “breaking” from the
diff’s shape alone (e.g. “this `pub fn`’s return type changed”). Instead,
trace the item’s visibility all the way up to the crate root:

- A `pub fn` on a `pub` type still isn’t part of the public API if it
  lives inside a module that isn’t `pub mod` the whole way up. For
  example, `s2::random` is declared as a private `mod random;`
  ([`src/s2/mod.rs`](src/s2/mod.rs)), so nothing inside it — no matter how
  it’s written — is reachable by a downstream crate; changing it is
  invisible to them.
- Conversely, changing something from `pub` to `pub(crate)` *removes* it
  from the public API even if nothing about its signature otherwise
  changed — that’s a real break, not an internal detail.

When describing a PR with both kinds of change, it’s worth splitting them
out explicitly (“breaking” vs. “purely internal, no external effect”)
rather than describing the diff as one flat list — it’s the difference
between something a downstream maintainer needs to act on and something
they’ll never notice.

## Release process 📦

Currently manual:

1. Open a PR bumping `version` in `Cargo.toml` per the SemVer rule above,
   informed by the categorized diff since the last tag (e.g. via `gh api
   repos/yjh0502/rust-s2/releases/generate-notes` against the previous tag).
2. Once merged, tag `master` as `vX.Y.Z` and push the tag.
3. The tag push triggers `.github/workflows/release.yml`, which runs the
   test suite, publishes to crates.io (via trusted OIDC publishing — no
   token to manage), and creates the GitHub Release with generated,
   categorized notes.
