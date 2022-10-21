# Prepare packages ----------------------------------------------------------------------------------------------------
import numpy as np
# ---------------------------------------------------------------------------------------------------------------------


# RANDOM FUNCTIONS ---------------------------------------------------------------------------------------------------

def normal_dist(mean, sd, sample_size, dataset_id=1):
    """
    @Definition:
                Each cross-section area point will have a normal distribution. A random value will be chosen from
                this distribution for each simulation. For example: for each of 5 simulations, 5 random values of
                chosen para will be picked from this normal distribution.

                This function is to generate a random normal distribution
    @References:
                None.
    @Arguments:
                mean (float):
                            This will be estimated value from LiDAR developed by Rose's codes
                sd (float):
                            Random value (but should be reasonable)
                sample_size (float):
                            Number of random values being generated from normal distribution
                dataset_id (int):
                            Ordinal number/ ID of random normal distribution of a chosen para
    @Returns:
                A random normal distribution
    """
    # Get seed
    np.random.seed(dataset_id)
    return np.random.normal(
        loc=mean,
        scale=sd,
        size=sample_size
    )

def normal_dist_number(
    mean, sd,
    sample_size, dataset_id_sample, dataset_id_number
):
    """
    @Definition:
                A function to generate a random number from a normal distribution
    @References:
                None.
    @Arguments:
                mean (float):
                            This will be estimated value from LiDAR developed by Rose's codes
                sd (float):
                            Random value (but should be reasonable)
                sample_size (float):
                            Number of random values being generated from normal distribution
                dataset_id_sample (int):
                            Ordinal number/ ID of random normal distribution of a chosen para
                dataset_id_number (int):
                            Ordinal number/ ID of random number chosen from random normal distribution
    @Returns:
                Random number from normal distribution
    """
    # Get array of normal distribution
    random_normal_arr = normal_dist(mean, sd, sample_size, dataset_id_sample)

    # Get seed
    np.random.seed(dataset_id_number)
    number = np.random.choice(random_normal_arr, replace=False)
    if number <= 0:
        number = number - np.min(random_normal_arr)
        if number == 0:
            number = 0.1
        else:
            pass
    else:
        pass
    return number

def uniform_dist_number(
    low, high, dataset_id
):
    """
    @Definition:
                A function to generate a random number from a uniform distribution
    @References:
                None.
    @Arguments:
                low an high (float):
                                Lower and upper limits of uniform distribution
                dataset_id (int):
                                ID of uniform distribution
    @Returns:
                A random number from uniform distribution
    """
    # Set up seed
    rng = np.random.default_rng(dataset_id)
    # Get random number
    number = rng.uniform(low, high, 1)[0]

    return number
# END RANDOM FUNCTIONS -----------------------------------------------------------------------------------------------
