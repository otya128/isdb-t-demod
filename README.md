# isdb-t-demod

2023/12/31のコミックマーケット103で頒布する合同誌Ariake Tea Partyに寄稿するにあたって実装したISDB-Tワンセグ復調の実装

詳細: <https://www.otyakai.xyz/c103/>

## 必要なもの

* 8.126984 MHzのサンプリングレートでI/Q信号を出力できるハードウェア (DTV02A-1T-Uなど)

## 使い方

ビルド:
```
cd isdb-t-demo
cargo build --release
```

DTV02A-1T-Uやそれに類するMSi001+MSi2500を使ったデバイスはデフォルトではVideo for Linux用のドライバが使われ、機能に制約があるためカーネルモジュールをアンロードするか`/etc/modprobe.d/blacklist.conf` に下記を追加するなどして使われないようにしたのちに`libusb`を用いる<https://github.com/f4exb/libmirisdr-4> もしくは <https://github.com/ericek111/libmirisdr-5> を使うことを推奨
```
blacklist msi001
blacklist msi2500
```

DTV02A-1T-Uを使い、473.142857 MHz (13ch)をサンプリング (帯域幅6MHz、zero-IF、12-bit):
```
miri_sdr -f 473142857 -s 8126984 -w 6000000 -i 0 -m 384_S16 > sample.raw
```

復調し、MPEG TS `test.ts`を出力:
```
./target/release/isdb-t-demod --input sample.raw --output test.ts
```

うまく復調できると以下のような出力が得られます。
```
CFO: -1
FRAME NOT DETECTED
...
FRAME DETECTED 0
symbol buffer not filled
FRAME DETECTED 1
symbol buffer not filled
FRAME DETECTED 2
CNR = 19.46 dB
byte buffer not filled
...
FRAME DETECTED 3
FRAME DETECTED 4
FRAME DETECTED 5
...
```

再生:
```
ffplay test.ts
```
