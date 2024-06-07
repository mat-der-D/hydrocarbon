# 概要
単結合のみを含む炭化水素の構造を列挙するプログラムである。
以下は main.rs 内の主な登場人物である。

- `do_search<N>(num_threads: usize, save_results: bool)`<br/>
  炭素数 `N` の炭化水素の構造を全列挙し, その数を標準出力する。<br/>
  `num_threads` は構造の同型判定の際に行われる並列計算のスレッド数を表す。<br/>
  `save_results` が true のときは, 各構造を表せる隣接行列をファイルに出力する。
  例えば `N` が 3 のとき, 以下の内容の .json ファイルが出力される。
  ```json
  [
    {
      "dim": 3,
      "matrix": "011101110",
      "hash": [
        0,
        0,
        0
      ]
    },
    {
      "dim": 3,
      "matrix": "011100100",
      "hash": [
        0,
        1,
        1
      ]
    }
  ]
  ```
  dim は行列の次元, matrix は行列の要素, hash は Morgan 法に基づいて計算されたハッシュ値である。
  matrix の部分はそれぞれ
  ```
  011
  101
  110
  ```
  と
  ```
  011
  100
  100
  ```
  を隣接行列とする構造があることを表す。すなわち各行列は長さ `N^2` のビット列として表現されている。
- `do_search_many!(..)`<br/>
  引数に与えた炭素数の探索を行うことを表す。例えば以下は炭素数 2~10 の探索を順に行う命令である。
  ```rust
  do_search_many!(2, 3, 4, 5, 6, 7, 8, 9, 10);
  ```
- `NUM_THREADS: usize`, `SAVE_RESULTS: bool`<br/>
  `do_search_many!(..)` の振る舞いの設定用の定数。それぞれ `do_search` の引数に与えられる。