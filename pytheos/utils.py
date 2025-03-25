# general utilies


def check_with_user() -> None:
    """
    Generic while loop to check if user wants to continue. Usually only called
    interally when something wrong has been detected within a method.
    """

    while True:
        answer = input("\nDo you want to continue? (yes/no): ").lower()
        if answer in ["y", "yes"]:
            print("Continuing...\n")
            break

        elif answer in ["n", "no"]:
            print("Aborting...\n")
            exit()
            break

        else:
            print("Invalid input. Please enter 'yes' or 'no'")

    return None
