yolo-fix:
    -__CARGO_FIX_YOLO=1 cargo fix --allow-dirty --all-features --all-targets --broken-code
    -__CARGO_FIX_YOLO=1 cargo clippy --fix --allow-dirty --all-features --all-targets --broken-code
    -cargo fmt