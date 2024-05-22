import click

VERSION = "0.1.0"

@click.command()
@click.version_option(VERSION, message="%(version)s")
@click.argument("paths", nargs=3)

def main(paths):
    print("Xenium main folder path: "+paths[0])
    print("Coordinates folder path: "+paths[1])
    print("Target folder path: "+paths[2])

if __name__ == "__main__":
    main()