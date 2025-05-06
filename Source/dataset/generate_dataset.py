import random

def generate_dataset(num_rows, num_nodes, weight_range=(1, 100)):
    with open("data.txt", "w") as f:
        for _ in range(num_rows):
            src = random.randint(0, num_nodes - 1)
            dest = random.randint(0, num_nodes - 1)
            while dest == src:
                dest = random.randint(0, num_nodes - 1)
            weights = [random.randint(*weight_range) for _ in range(5)]
            row = [src, dest] + weights
            f.write(','.join(map(str, row)) + ',\n')  # trailing comma and newline
        f.write("~")

# Example usage
generate_dataset(num_rows=10000000, num_nodes=10000)
