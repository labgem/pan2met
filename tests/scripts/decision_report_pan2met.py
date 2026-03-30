import sys


def main():
    with open(sys.argv[1], "r") as input_file:
        with open(sys.argv[2], "w") as output_file:
            for row in input_file:
                if row.startswith("DECISION:"):
                    pathway = row.replace("DECISION:", "").replace("\n", "")
                elif row.startswith("REJECT:"):
                    reason = row.replace("REJECT:", "").replace("\n", "")
                    decision = False
                elif row.startswith("ACCEPT:"):
                    reason = row.replace("ACCEPT:", "").replace("\n", "")
                    decision = True
                if row.startswith("REJECT") or row.startswith("ACCEPT"):
                    output_file.write(
                        "\t".join([pathway, str(decision), reason]) + "\n"
                    )


if __name__ == "__main__":
    main()
