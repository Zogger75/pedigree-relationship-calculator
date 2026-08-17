Pedigree Relationship Calculator

This project demonstrates recursive pedigree traversal, coefficient of parentage calculation, relationship matrix generation, and exportable reporting for plant breeding datasets.

COP calculations are based on the methodology described by Graham McLaren, Ian DeLacy, and Jose Crossa in Routine Computation and Visualization of Coefficients of Parentage Using the International Crop Information System (ICIS).

## Input Pedigree File Format

The tool reads a pipe-delimited text file with one germplasm record per row.

Required columns:

| Column | Description |
|---|---|
| `gid` | Unique numeric identifier for the germplasm line |
| `name` | Human-readable germplasm or line name |
| `dam` | Female parent / first parent GID. Use `0` if unknown |
| `sire` | Male parent / second parent GID. Use `0` if unknown |
| `method_code` | Breeding method used to derive the line |
| `include` | Use `Y` to include the line in the final COP matrix; use `N` for ancestors only |

### Method Codes

| Code | Meaning | How the calculator treats it |
|---|---|---|
| `FND` | Founder line with no known parents | Treated as a base ancestor |
| `GEN` | Generative cross | Uses both `dam` and `sire` as genetic parents |
| `DER` | Derived or selected line | Uses the listed parentage |
| `SEL` | Selection from an existing line or population | Treated similarly to `DER` |
| `MAN` | Management, increase, release, or naming step | Treated as no genetic change from the source parent |

Example:

```text
gid|name|dam|sire|method_code|include
1001|Founder_A|0|0|FND|N
2001|Cross_AB|1001|1002|GEN|N
3001|Selected_Line_1|2001|0|DER|Y
