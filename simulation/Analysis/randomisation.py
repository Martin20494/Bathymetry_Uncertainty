# Prepare packages  ----------------------------------------------------------------------------------------------------
import numpy as np                                        # For array manipulation
# ----------------------------------------------------------------------------------------------------------------------
def random_values_generation(num_seed=1,
                             num_simulation=300,
                             rotation_boundaries=None,
                             translation_boundaries_x=None,
                             translation_boundaries_y=None,
                             type_random='uniform',
                             replacement=False):
    """
    @Definition:
                A function to generate random numbers for angle, x, and y. If 'systematic', the values for rotation
                and translation will be chosen with interval. If 'uniform', the values for rotation and translation
                will be chosen randomly without replacement.
    @References:
                https://florimond.dev/en/posts/2018/08/python-mutable-defaults-are-the-source-of-all-evil/
                https://stackoverflow.com/questions/66966057/setting-a-maximum-number-of-replacements-for-sample-with-numpy-choice
                https://numpy.org/doc/stable/reference/random/generated/numpy.random.Generator.choice.html#numpy.random.Generator.choice
                https://numpy.org/doc/stable/reference/random/generated/numpy.random.Generator.integers.html#numpy.random.Generator.integers
                https://numpy.org/doc/stable/reference/random/generator.html
    @Arguments:
                num_seed (int):
                                            The ordinal number of seed that identifies specific data. Default is 1
                num_simulation (int):
                                            The quantity of simulations. Default is 300
                rotation_boundaries (list):
                                            A list of values representing the boundaries of random angles
                                            If 'systematic', for instance, 'rotation_boundaries = [0, 91, 5]' means the
                                            angle will be chosen between 0 and 90 with interval of 5 degree.
                                            If 'uniform', for instance, 'rotation_boundaries = [0, 91]' means the angle
                                            will be randomly chosen between 0 and 90.
                                            Default is [0, 91, 5]
                translation_boundaries_x (list):
                                            A list with two values representing the boundaries of random x.
                                            If 'systematic', for instance, 'translation_boundaries_x = [0, 6, 1]' means
                                            x will be chosen between 0 and 5 meters with interval of 1 meter.
                                            If 'uniform', for instance, 'translation_boundaries_x = [0, 6]' means x will
                                            be randomly chosen between 0 and 5 meters.
                                            Default is [0, 6, 1]
                translation_boundaries_y (list):
                                            A list with two values representing the boundaries of random y.
                                            If 'systematic', for instance, 'translation_boundaries_y = [0, 6, 1]' means
                                            y will be chosen between 0 and 5 meters with interval of 1 meter.
                                            If 'uniform', for instance, 'translation_boundaries_y = [0, 6]' means y will
                                            be chosen between 0 and 6 meters.
                                            Default is [0, 6, 1]
                type_random (string):
                                            Could be 'systematic' or 'uniform'
                replacement (boolean):
                                            Could be False or True to identify uniform randomisation
                                            without/with replacement
    @Return:
                array_transformation (array):
                                            An array of random angle, x, and y
    -----------
    """
    # Set defaults for boundaries of rotation and translations ---------------------------------------------------------
    if rotation_boundaries is None:
        rotation_boundaries = [0, 91, 5]
    elif translation_boundaries_x is None:
        translation_boundaries_x = [0, 6, 1]
    elif translation_boundaries_y is None:
        translation_boundaries_y = [0, 6, 1]
    # ------------------------------------------------------------------------------------------------------------------

    # Set up 'systematic' randomisation --------------------------------------------------------------------------------
    if type_random == 'systematic':
        # Get empty list
        list_system = []

        # Generate random angle, x, and y under systematic randomization
        for angle_system in range(rotation_boundaries[0], rotation_boundaries[1], rotation_boundaries[2]):
            for x_system in range(translation_boundaries_x[0], translation_boundaries_x[1],
                                  translation_boundaries_x[2]):
                for y_system in range(translation_boundaries_y[0], translation_boundaries_y[1],
                                      translation_boundaries_y[2]):
                    new_system = [angle_system, x_system, y_system]
                    list_system.append(new_system)

        # Convert list into array
        array_transformation = np.array(list_system)
    # ------------------------------------------------------------------------------------------------------------------

    # Set up 'uniform' randomisation -----------------------------------------------------------------------------------
    else:
        # Set up seed
        rng = np.random.default_rng(num_seed)

        # Get empty list
        list_uni = []

        # Generate random angle, x, and y under uniform randomization
        for i in range(num_simulation):
            while True:
                angle_uni = rng.uniform(low=rotation_boundaries[0], high=rotation_boundaries[1], size=1)[0]
                x_uni = rng.uniform(low=translation_boundaries_x[0], high=translation_boundaries_x[1], size=1)[0]
                y_uni = rng.uniform(low=translation_boundaries_y[0], high=translation_boundaries_y[1], size=1)[0]
                new_uni = [angle_uni, x_uni, y_uni]

                # Set up 'not replacement' if replacement == False
                if replacement is True:
                    list_uni.append(new_uni)
                    break
                else:
                    # Eliminate the duplicates
                    if new_uni in list_uni:
                        continue
                    else:
                        list_uni.append(new_uni)
                    break
        # --------------------------------------------------------------------------------------------------------------

        # Convert list into array
        array_transformation = np.array(list_uni)

    return array_transformation