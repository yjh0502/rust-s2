# rust-s2 🌐

[![test coverage](https://coveralls.io/repos/github/yjh0502/rust-s2/badge.svg?branch=master)](https://coveralls.io/github/yjh0502/rust-s2?branch=master)
[![crates.io](https://img.shields.io/crates/v/s2)](https://crates.io/crates/s2)
[![docs](https://docs.rs/s2/badge.svg)](https://docs.rs/s2/latest/s2/)
[![MSRV](https://img.shields.io/crates/msrv/s2)](https://crates.io/crates/s2)

A Rust implementation of [Google’s S2 geometry
library](https://github.com/google/s2geometry): spherical geometry on the
surface of a sphere, built around a hierarchical decomposition of the
sphere into cells. It’s what you reach for when flat-plane geometry isn’t
good enough — proximity search, region containment, and indexing for
points, lines, and polygons on the globe.

## Status 🚧

This is an ongoing port, further along in some areas than others — in
particular, lines and polygons aren’t implemented yet. Rather than
maintain a status table here that quietly goes stale, have a look at the
[issue tracker](https://github.com/yjh0502/rust-s2/issues) for what’s
implemented, what’s missing, and what’s being worked on.

## Contributing 🙌

Contributions of any size are very welcome — typo fixes, bug reports,
new tests, or porting a whole piece of functionality. See
[CONTRIBUTING.md](CONTRIBUTING.md) to get started; issues labeled
[`good first
issue`](https://github.com/yjh0502/rust-s2/labels/good%20first%20issue)
are a good place to look if you’re new here.

Please also take a look at our [Code of Conduct](CODE_OF_CONDUCT.md), and
our [Security Policy](SECURITY.md) if you need to report a vulnerability.

## Credits ❤️

This library exists thanks to [@yjh0502](https://github.com/yjh0502), who
started and has driven most of this port, with substantial contributions
from [@brawer](https://github.com/brawer), [@nside](https://github.com/nside),
[@Louisvranderick](https://github.com/Louisvranderick), and
[@HLennart](https://github.com/HLennart) — and everyone else in the
[contributors graph](https://github.com/yjh0502/rust-s2/graphs/contributors).
Thank you all! 🎉

## License ⚖️

Apache-2.0 — see [LICENSE](LICENSE).
