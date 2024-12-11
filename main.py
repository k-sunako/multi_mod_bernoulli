from src.bernoulli import rational
import sys
sys.set_int_max_str_digits(0)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('n', type=int)
    parser.add_argument('--to_file', action='store_true')

    args = parser.parse_args()
    Bn = rational(args.n)
    if args.to_file:
        with open('out.txt', 'w') as f:
            f.write(str(Bn.numerator))
            f.write('\n')
            f.write(str(Bn.denominator))
    else:
        print(Bn)
