# general utilies


def check_with_user() -> None:
    """
    Generic while loop to check if user wants to continue. Usually only called
    interally when something wrong has been detected within a method.
    """

    while True:
        answer = input("Do you want to continue? (yes/no): ").lower()
        if answer in ["y", "yes"]:
            print("Continuing...")
            break

        elif answer in ["n", "no"]:
            print("Aborting...")
            exit()
            break

        else:
            print("Invalid input. Please enter 'yes' or 'no'")

    return None


def sum(numbers=list) -> float:
    # only for learning how to write unit tests...
    total = 0
    for number in numbers:
        total += number

    return total
