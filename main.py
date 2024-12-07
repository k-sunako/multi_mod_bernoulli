from src.bernoulli import rational

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('n', type=int)

    args = parser.parse_args()
    print(rational(args.n))
