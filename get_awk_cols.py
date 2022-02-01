import sys

def get_col(arg):
    total = 0
    for char in arg:
        total *= 26
        total += ord(char) - ord('A') + 1
    return total


def get_cols(arg):
    arg = arg.upper()
    arg = arg.replace('_', '-')
    bounds = arg.split('-')
    output = ''
    if len(bounds) == 1:
        output = f'${get_col(bounds[0])}'
    elif len(bounds) == 2:
        for col in range(get_col(bounds[0]), get_col(bounds[1]) + 1):
            output += f'${col},'
    else:
        assert False
    return output.strip(',')


def get_output(arg_list):
    if type(arg_list) == type('string'):
        arg_list = arg_list.split()
    output = ''
    for arg in arg_list:
        output += get_cols(arg) + ','
    return output.strip(',')


def main():
    print(get_output(sys.argv[1:]))


if __name__ == '__main__':
    main()
