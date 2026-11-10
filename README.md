## LiftOver of Repli-seq WaveSignal tracks from hg19 to hg38

This workflow converts 16 UCSC ENCODE Repli-seq WaveSignal tracks from the hg19 human genome assembly to hg38. The commands are for an Apple Silicon Mac with Homebrew installed. Run them in order from the same Terminal window.

### Step 1 — Install the required tools

```bash
brew install wget parallel xz bedtools
```

### Step 2 — Create a directory for the UCSC executables

```bash
mkdir -p "$HOME/bin"
```

### Step 3 — Download the UCSC executables for Apple Silicon

```bash
wget -O "$HOME/bin/liftOver" \
    https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/liftOver

wget -O "$HOME/bin/bigWigToBedGraph" \
    https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/bigWigToBedGraph

wget -O "$HOME/bin/bedGraphToBigWig" \
    https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/bedGraphToBigWig
```

### Step 4 — Make the UCSC tools executable

```bash
chmod +x "$HOME/bin/liftOver" \
    "$HOME/bin/bigWigToBedGraph" \
    "$HOME/bin/bedGraphToBigWig"
```

### Step 5 — Add the UCSC tools directory to your PATH

Run this once if `$HOME/bin` is not already configured in your `~/.zshrc`.

```bash
echo 'export PATH="$HOME/bin:$PATH"' >> "$HOME/.zshrc"
```

### Step 6 — Load the configuration and use Bash for Parallel jobs

```bash
source "$HOME/.zshrc"
export PARALLEL_SHELL=/bin/bash
```

### Step 7 — Set up the working folder and download the reference files

The chain file maps hg19 coordinates to hg38. The chromosome sizes file is used to build the final hg38 bigWigs.

```bash
mkdir -p "$HOME/res/bin"
cd "$HOME/res/bin"

curl -fL -o hg19ToHg38.over.chain.gz \
    https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

gzip -t hg19ToHg38.over.chain.gz

curl -fL -o hg38.chrom.sizes \
    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
```

### Step 8 — Download and count the hg19 WaveSignal tracks

```bash
wget -r -l1 -np -nd -c \
    -A '*WaveSignal*.bigWig' \
    -P repli_seq_processed \
    https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/

find repli_seq_processed -type f -name '*WaveSignal*.bigWig' -size +0c | wc -l
```

The expected input count is **16**.

### Step 9 — Create the output directories

```bash
mkdir -p outbed outbw unmapped
```

### Step 10 — Convert bigWig to bedGraph and map coordinates to hg38

Mapped intervals are saved in `outbed/`. Intervals that liftOver could not map are recorded in `unmapped/`.

```bash
find repli_seq_processed -type f -name '*WaveSignal*.bigWig' -size +0c -print0 |
    parallel -0 --jobs 4 '
        set -o pipefail
        liftOver <(bigWigToBedGraph {} /dev/stdout) \
            hg19ToHg38.over.chain.gz \
            /dev/stdout unmapped/{/.}_unmapped.txt -bedPlus=4 |
            LC_ALL=C sort -k1,1 -k2,2n > outbed/{/.}.bedGraph
    '
```

### Step 11 — Create temporary interval templates

For each track, collect the start and end coordinates, sort them, remove duplicates, and use consecutive boundaries to define smaller intervals.

```bash
find outbed -type f -name '*WaveSignal*.bedGraph' \
    ! -name '*.hg38.sorted.bedGraph' -size +0c -print0 |
    parallel -0 --jobs 4 '
        set -o pipefail
        (cut -f2 {}; cut -f3 {}) |
            sort -nu |
            perl -ne '\''BEGIN{$i=0} chomp; push @F, $_; if($i++){print "chr_placeholder\t$F[$i-2]\t$F[$i-1]\n"}'\'' \
            > {.}_master.tmp
    '
```

### Step 12 — Split tracks by chromosome and resolve overlapping intervals

Process chromosomes 1–22, X and Y. When multiple mapped intervals cover the same segment, keep the highest signal value.

```bash
find outbed -type f -name '*WaveSignal*_master.tmp' -size +0c -print0 |
    parallel -0 --jobs 4 '
        set -o pipefail
        TEMPLATE_PATH={}
        DATA_FILE="${TEMPLATE_PATH%_master.tmp}.bedGraph"
        TEMPLATE_NAME=${TEMPLATE_PATH##*/}

        for chr in {1..22} X Y
        do
            TILE_FILE="outbed/${TEMPLATE_NAME%.tmp}_b.${chr}.bed"

            sed "s/chr_placeholder/chr${chr}/" "$TEMPLATE_PATH" \
                > "$TILE_FILE" || exit 1

            bedtools intersect -a "$DATA_FILE" -b "$TILE_FILE" |
                LC_ALL=C sort -k2,2n -k3,3n -k4,4nr |
                perl -lane '\''print unless $h{$F[0],$F[1],$F[2]}++'\'' \
                > "${DATA_FILE}..${chr}" || exit 1
        done
    '
```

### Step 13 — Combine and sort the chromosome files for each track

```bash
find outbed -type f -name '*WaveSignal*.bedGraph' \
    ! -name '*.hg38.sorted.bedGraph' -size +0c -print0 |
    parallel -0 --jobs 4 '
        set -o pipefail
        input={}
        output="${input%.bedGraph}.hg38.sorted.bedGraph"

        cat "$input"..{1..22} "$input"..X "$input"..Y |
            LC_ALL=C sort -k1,1 -k2,2n -k3,3n > "$output"
    '
```

### Step 14 — Convert the sorted hg38 bedGraph files to bigWig

```bash
find outbed -type f -name '*.hg38.sorted.bedGraph' -size +0c -print0 |
    parallel -0 --jobs 4 '
        input={}
        name="${input##*/}"

        bedGraphToBigWig "$input" hg38.chrom.sizes \
            "outbw/${name%.hg38.sorted.bedGraph}.hg38.bigWig"
    '
```

### Step 15 — List and count the final hg38 bigWig files

```bash
find outbw -type f -name '*.hg38.bigWig' -size +0c | sort
find outbw -type f -name '*.hg38.bigWig' -size +0c | wc -l
```

The expected result is **16 non-empty hg38 bigWig files** in `outbw/`. The `outbed/` directory contains intermediate files.
