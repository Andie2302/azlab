#!/usr/bin/env bash
# Prüfmatrix für den azlab-Workspace (siehe todo/05-qualitaet.md).
#
# Schreibt das Ergebnis nach pruefung/ergebnis.txt.
#
# Aufruf:  scripts/pruefen.sh
# Dauer:   erster Lauf voraussichtlich wenige Minuten, danach deutlich kürzer.

set -u
cd "$(dirname "$0")/.." || exit 1

# Ohne inkrementellen Compiler-Cache prüfen: Ein beschädigter Cache hat schon
# einmal einen internen Compilerfehler ausgelöst, der nichts mit dem Code zu
# tun hatte. Für die kleine Bibliothek kostet das nur Sekunden.
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

for combo in "${combos[@]}"; do
    # Absichtlich ohne Anführungszeichen: $combo enthält mehrere Argumente.
    # shellcheck disable=SC2086
    run cargo test -p azmath $combo
    # shellcheck disable=SC2086
    run cargo clippy -p azmath --all-targets $combo -- -D warnings
done

run cargo doc -p azmath --no-deps --all-features

# Echter no_std-Beweis: Build für ein Ziel ohne Betriebssystem.
target="thumbv7em-none-eabihf"
if rustup target list --installed 2> /dev/null | grep -qx "$target"; then
    run cargo build -p azmath --target "$target" --no-default-features
    run cargo build -p azmath --target "$target"
else
    echo "== $target nicht installiert, übersprungen (rustup target add $target)" | tee -a "$out"
fi

echo | tee -a "$out"
if [ "$fail" -eq 0 ]; then
    echo "ERGEBNIS: alles OK" | tee -a "$out"
else
    echo "ERGEBNIS: FEHLER, Details in $out" | tee -a "$out"
fi
exit "$fail"
