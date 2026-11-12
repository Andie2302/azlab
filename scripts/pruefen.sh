#!/usr/bin/env bash
# Prüfmatrix für den azlab-Workspace (siehe todo/05-qualitaet.md).
#
# Schreibt das Ergebnis nach pruefung/ergebnis.txt.
#
# Aufruf:  scripts/pruefen.sh
# Dauer:   erster Lauf voraussichtlich einige Minuten (die Anzeige für die
#          Beispiele bringt größere Abhängigkeiten mit), danach deutlich kürzer.

set -u
cd "$(dirname "$0")/.." || exit 1

# Ohne inkrementellen Compiler-Cache prüfen: Ein beschädigter Cache hat schon
# einmal einen internen Compilerfehler ausgelöst, der nichts mit dem Code zu
# tun hatte.
export CARGO_INCREMENTAL=0

out="pruefung/ergebnis.txt"
mkdir -p pruefung
: > "$out"
fail=0

run() {
    echo "== $*" | tee -a "$out"
    if "$@" >> "$out" 2>&1; then
        echo "   OK" | tee -a "$out"
    else
        echo "   FEHLER" | tee -a "$out"
        fail=1
    fi
}

{
    date
    rustc --version
    cargo --version
    echo
} >> "$out" 2>&1

# Jede Feature-Kombination muss bauen, Tests bestehen und Clippy-frei sein.
combos=(
    "--no-default-features"
    "--no-default-features --features libm"
    "--no-default-features --features std"
    "--no-default-features --features alloc"
    "--all-features"
)

for crate in azmath azphysics; do
    for combo in "${combos[@]}"; do
        # Absichtlich ohne Anführungszeichen: $combo enthält mehrere Argumente.
        # shellcheck disable=SC2086
        run cargo test -p "$crate" $combo
        # shellcheck disable=SC2086
        run cargo clippy -p "$crate" --all-targets $combo -- -D warnings
    done
done

# Die Anzeige für die Beispiele braucht std und hat keine Feature-Schalter.
run cargo test -p azviz
run cargo clippy -p azviz --all-targets -- -D warnings

run cargo doc --workspace --no-deps --all-features

# Echter no_std-Beweis: Build für ein Ziel ohne Betriebssystem.
target="thumbv7em-none-eabihf"
if rustup target list --installed 2> /dev/null | grep -qx "$target"; then
    for crate in azmath azphysics; do
        run cargo build -p "$crate" --target "$target" --no-default-features
        run cargo build -p "$crate" --target "$target"
    done
    # Simulationen mit Heap, aber ohne Betriebssystem.
    run cargo build -p azphysics --target "$target" --features alloc
else
    echo "== $target nicht installiert, übersprungen (rustup target add $target)" | tee -a "$out"
fi

# Die zugesagte Mindestversion (rust-version im Workspace-Manifest): Die
# beiden Bibliotheken müssen sich mit ihr übersetzen lassen. Die Anzeige für
# die Beispiele ist davon ausgenommen; ihre Abhängigkeiten verlangen mehr.
msrv="$(sed -n 's/^rust-version *= *"\(.*\)"/\1/p' Cargo.toml | head -n 1)"
if [ -n "$msrv" ] && rustup toolchain list 2> /dev/null | grep -q "^$msrv"; then
    run cargo "+$msrv" check -p azmath -p azphysics --no-default-features
    run cargo "+$msrv" check -p azmath -p azphysics --all-features
else
    echo "== Rust $msrv nicht installiert, übersprungen (rustup toolchain install $msrv --profile minimal)" | tee -a "$out"
fi

echo | tee -a "$out"
if [ "$fail" -eq 0 ]; then
    echo "ERGEBNIS: alles OK" | tee -a "$out"
else
    echo "ERGEBNIS: FEHLER, Details in $out" | tee -a "$out"
fi
exit "$fail"
