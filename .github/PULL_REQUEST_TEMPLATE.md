<!-- Thanks for the PR! A couple of things to double check before/while it's open — see CONTRIBUTING.md for the full detail on any of these. -->

- [ ] PR title follows [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/) — it becomes the commit message on `master`, since PRs are squash-merged.
- [ ] `cargo fmt --all`, `cargo clippy --all-features --all-targets`, and `cargo test --all-features` pass locally.
- [ ] If this ports or fixes S2 algorithm code: checked against the [C++ source](https://github.com/google/s2geometry) (not just Go), and ported its tests too.

<!-- Delete anything above that doesn't apply. -->
