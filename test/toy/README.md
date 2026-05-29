## GapFill notes

To convert data
- `toy2dat.py`  converts Amy's toy problem .xlsx data into a 'known flux' format file (`toyflux.dat`)

- `tsvcvt.py` takes reaction and compound .tsv (Tab Separated Values) files, and produces `all.dat` - a list of all reactions and compounds

- `makedat` takes known flux file `toyflux.dat` and all reactions file `all.dat` and produces `toy.dat`, the input file for `gapfill`

  

To run gapfill with known flux

```bash
gapfill -s 1 -v CTS	-f toyflux.dat -o sol.out -g sol.dig toy.dat
```

Notes:

- -s is seed - can be other numbers
- CTS is the continuous solver
- Known flux in `toyflux.dat`
- Outputs `sol.out` and `sol.dig`



