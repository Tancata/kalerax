#/bin/sh

ask_yes_no() {
    while true; do
        read -p "$1 (y/N): " answer
        case $answer in
            [Yy]* ) return 0 ;;
            ""|[Nn]* ) return 1 ;;
            * ) echo "Please answer yes or no." ;;
        esac
    done
}

num_cores=$(nproc 2>/dev/null || echo 1)

# Apply the kalerax patches to the GeneRaxCore submodule (idempotent).
for p in patches/generaxcore_consensus_export_guards.patch \
         patches/generaxcore_wgd_export_synteny.patch; do
  [ -f "$p" ] || continue
  if (cd ext/GeneRaxCore && patch -p1 -N --dry-run --silent < "../../$p" > /dev/null 2>&1); then
    (cd ext/GeneRaxCore && patch -p1 -N --silent < "../../$p") && echo "applied $p"
  else
    echo "skipping $p (already applied or not applicable)"
  fi
done

# Example usage:
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release

if ask_yes_no "Build only kalerax binary (otherwise build everything)?"; then
    echo "Building kalerax..."
    make -j $num_cores kalerax
else
    echo "Building everything..."
    make -j $num_cores
fi
