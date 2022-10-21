# Prepare packages -----------------------------------------------------------------------------------------------------
import numpy as np                                                   # For data manipulation
import matplotlib.colors as mcolors                                  # For coloring
# ----------------------------------------------------------------------------------------------------------------------

def hex_to_rgb(value_func):
    """
    @Definition:
                A function to convert hex to red green blue (rgb) colors. The code is based on the link in reference
                (by Kerry Halupka)
    @References:
                https://towardsdatascience.com/beautiful-custom-colormaps-with-matplotlib-5bab3d1f0e72
    @Arguments:
                value_func (string):
                                 Hex color code under 6-characters string format
    @Returns:
                rgb_value_func (tuple):
                                 A tuple of rgb values with length of 3
    """
    # Remove '#' in string
    hex_value_func = value_func.strip('#')

    # Get the quantity of colors
    level_func = len(hex_value_func)

    # Convert color from hex to rgb and store in a tuple
    rgb_value_func = tuple(
        int(hex_value_func[i:i + level_func // 3], 16) for i in range(0, level_func, level_func // 3))

    return rgb_value_func

def rgb_to_dec(value_func):
    """
    @Definition:
                A function to convert rgb to decimal colors (dividing each value by 256). The code is based on the link
                in reference (by Kerry Halupka)
    @References:
                https://towardsdatascience.com/beautiful-custom-colormaps-with-matplotlib-5bab3d1f0e72
    @Arguments:
                value_func (tuple):
                                 A tuple of rgb color code (from 0 to 256) with length of 3
    @Returns:
                dec_value_func (list):
                                A list of color decimal values with length of 3
    """
    return [val / 256 for val in value_func]

def get_gradient_cmap(hex_list_func, float_list_func=None):
    """
    @Definition:
                A function to create gradient colors. The code is based on the link in reference (by Kerry Halupka)
                If float_list_func is None, colour map will graduate linearly between each color in hex_list
                If float_list_func is not None, each color in hex_list_func is mapped to the respective location in float_list_func
    @References:
                https://towardsdatascience.com/beautiful-custom-colormaps-with-matplotlib-5bab3d1f0e72
    @Arguments:
                hex_list_func (list):
                                                A list of hex code color under string format
                float_list_func (list):
                                                A list of floats (between 0 and 1), same length as hex_list_func.
                                                Must start with 0 and end with 1.
    @Returns:
                colour_map (color map in matplotlib):
                                                Color under matplotlib color map
    """
    # Get rgb list
    rgb_list = [rgb_to_dec(hex_to_rgb(color_code)) for color_code in hex_list_func]

    # Check float list
    if float_list_func:
        pass
    else:
        float_list_func = list(np.linspace(0, 1, len(rgb_list)))

    # Build up gradient colors
    color_dict = dict()
    for number, color in enumerate(['red', 'green', 'blue']):
        color_list = [[float_list_func[i], rgb_list[i][number], rgb_list[i][number]] for i in
                      range(len(float_list_func))]
        color_dict[color] = color_list

    color_map = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=color_dict, N=256)

    return color_map


### LIST OF HEX COLORS #################################################################################################

hex_list0 = ['#0091ad', '#3fcdda', '#83f9f8', '#d6f6eb', '#fdf1d2', '#f8eaad', '#faaaae', '#ff57bb']
hex_list1 = ["#8fe2ff", "#8dbdff", "#5999ff", "#3d50ff", "#6b01ff", "#4800c4", "#e200ff"]
hex_list2 = ["#0466C8", "#0353A4", "#023E7D", "#002855", "#001845", "#001233", "#33415C", "#5C677D", "#7D8597",
             "#979DAC"]
hex_list3 = ["#B7094C", "#A01A58", "#892B64", "#723C70", "#5C4D7D", "#455E89", "#2E6F95", "#1780A1", "#0091AD"]
hex_list4 = ["#25CED1", "#FFFFFF", "#FCEADE", "#FF8A5B", "#EA526F"]
hex_list5 = ["#EA698B", "#D55D92", "#C05299", "#AC46A1", "#973AA8", "#822FAF", "#6D23B6", "#6411AD", "#571089",
             "#47126B"]
hex_list6 = ["#25CED1", "#00B2CA", "F7C59F", "#FF8A5B", "#086375"]
hex_list7 = ["#B76935", "#A56336", "#935E38", "#815839", "#6F523B", "#5C4D3C", "#4A473E", "#38413F", "#263C41",
             "#143642"]
hex_list8 = ["#522E38", "#E05780", "#FF9EBB", "#A06CD5", "#3C096C", "#B298DC", "#979DAC", "#5C677D", "#0D1B2A"]
hex_list9 = ["#FF595E", "#FFCA3A", "#8AC926", "#1982C4", "#6A4C93"]
hex_list10 = ["#390099", "#9E0059", "#FF0054", "#FF5400", "#FFBD00"]
hex_list11 = ["#70D6FF", "#FF70A6", "#FF9770", "#FFD670", "#E9FF70"]
hex_list12 = ["#001219", "#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03", "#AE2012",
              "#9B2226"]
hex_list13 = ["#004B23", "#006400", "#007200", "#008000", "#38B000", "#70E000", "#9EF01A", "#CCFF33"]
hex_list14 = ["#90F1EF", "#FFD6E0", "#FFEF9F", "#C1FBA4", "#7BF1A8"]
hex_list15 = ["#3D348B", "#7678ED", "#F7B801", "#F18701", "#F35B04"]
hex_list16 = ["#54478C", "#2C699A", "#048BA8", "#0DB39E", "#16DB93", "#83E377", "#B9E769", "#EFEA5A", "#F1C453",
              "#F29E4C"]
hex_list17 = ["#EE6055", "#60D394", "#AAF683", "#FFD97D", "#FF9B85"]
hex_list18 = ["#00A6FB", "#0582CA", "#006494", "#003554", "#051923"]
hex_list19 = ["#f79256ff", "#fbd1a2ff", "#7dcfb6ff", "#00b2caff", "#1d4e89ff"]
hex_list20 = ["#083d77ff", "#ebebd3ff", "#f4d35eff", "#ee964bff", "#f95738ff"]
hex_list21 = ["#1d4e89ff", "#fbd1a2ff", "#7dcfb6ff", "#00b2caff", "#1d4e89ff"]
hex_list22 = ["#072ac8ff", "#1e96fcff", "#a2d6f9ff", "#fcf300ff", "#ffc600ff"]
hex_list23 = ["#ff499eff", "#d264b6ff", "#a480cfff", "#779be7ff", "#49b6ffff"]
hex_list24 = ["#d62839ff", "#ba324fff", "#175676ff", "#4ba3c3ff", "#cce6f4ff"]
hex_list25 = ["2d728f", "3b8ea5", "f5ee9e", "f49e4c", "ab3428"]
hex_list26 = ["9b5de5", "f15bb5", "fee440", "00bbf9", "00f5d4"]
hex_list27 = ["ef476f", "ffd166", "06d6a0", "118ab2", "073b4c"]
hex_list28 = ["001219", "005f73", "0a9396", "94d2bd", "e9d8a6", "ee9b00", "ca6702", "bb3e03", "ae2012", "9b2226"]
hex_list29 = ["2d00f7", "6a00f4", "8900f2", "a100f2", "b100e8", "bc00dd", "d100d1", "db00b6", "e500a4", "f20089"]
hex_list30 = ["54478c", "2c699a", "048ba8", "0db39e", "16db93", "83e377", "b9e769", "efea5a", "f1c453", "f29e4c"]
hex_list31 = ["2d728f", "3b8ea5", "f5ee9e", "f49e4c", "ab3428"]
hex_list32 = ["d62839", "ba324f", "175676", "4ba3c3", "cce6f4"]
hex_list33 = ["ff499e", "d264b6", "a480cf", "779be7", "175676"]
hex_list34 = ["0c0a3e", "7b1e7a", "b33f62", "f9564f", "f3c677"]
hex_list35 = ["ee6055", "60d394", "aaf683", "ffd97d", "ff9b85"]
hex_list36 = ["247ba0", "70c1b3", "b2dbbf", "f3ffbd", "ff1654"]
hex_list37 = ["03071e", "370617", "6a040f", "9d0208", "d00000", "dc2f02", "e85d04", "f48c06", "faa307", "ffba08"]
hex_list38 = ["3d5a80", "98c1d9", "e0fbfc", "ee6c4d", "293241", "#e60000"]
hex_list39 = ["8ecae6", "219ebc", "023047", "ffb703", "fb8500", "#54D852"]
hex_list40 = ["6411ad", "f15bb5", "fee440", "F9C74F", "00bbf9", "184E77"]
hex_list41 = ["0c0a3e", "7b1e7a", "b33f62", "f9564f", "f3c677", "#e60000"]
hex_list42 = ["2d728f", "3b8ea5", "f5ee9e", "f49e4c", "ab3428", "#e60000"]
hex_list43 = ["6e44ff", "b892ff", "ffc2e2", "ff90b3", "BA324F", "#e60000"]
hex_list44 = ["000814", "001d3d", "003566", "ffc300", "ffd60a", "#e60000"]
hex_list45 = ["001427", "708d81", "f4d58d", "bf0603", "8d0801", "#e60000"]
hex_list46 = ["#e65154", "#26b6ff", "#67e6d1", "#cd76d6", "#ffca8c", "#fff2b3", "#ff8cd9", "#d99d5b", "#c8f2a9",
              "#d4b8ff"]
hex_list47 = ["#3900b3", "#714dbf", "#9e6b90", "#cf9270", "#ebb698"]
hex_list48 = ["#6690ff", "#526aad", "#423b38", "#945e4c", "#ff9573"]
hex_list49 = ["#0051f2", "#6da7fb", "#daffa5", "#d9963a", "#b35a00"]
hex_list50 = ["#65a5ff", "#4366c4", "#5b2e73", "#b98600", "#ffb900"]
hex_list51 = ["#23ccff", "#2c8eac", "#474333", "#9b8020", "#ffc800"]
hex_list52 = ["#00aaff", "#0072aa", "#1c4d31", "#bd7e00", "#ffaa00", "#990901"]
hex_list53 = ["#954151", "#db9793", "#ead98b", "#7295bf", "#1d4e89", "#e60000"]
hex_list54 = ["#0080ff", "#005bb2", "#472459", "#c40f0f", "#ff4d4d", "#e60000"]
hex_list55 = ["#a00000", "#df7b7b", "#dae695", "#47b0df", "#008fbf", "#e60000"]
hex_list56 = ["#a00000", "#df7b7b", "#dae695", "#47b0df", "#008fbf", "#e60000"]
hex_list57 = ["#ff7040", "#bc3f2b", "#592e73", "#1a99c0", "#33ddff"]
hex_list58 = ["#00c3ff", "#4e858c", "#47443a", "#94593b", "#ff4040"]
hex_list59 = ["#ffea8c", "#b3ab60", "#4b595e", "#6693c8", "#aadbff"]
hex_list60 = ["#ed5151", "#149ece", "#a7c636", "#9e559c", "#fc921f", "#ffde3e"]
hex_list61 = ["#c65a18", "#f36f20", "#f7975e", "#fbc09b", "#fdd4ba"]
hex_list62 = ["#c65a18", "#f7975e", "#fdd4ba"]
hex_list63 = ["#ffff00", "#a19700", "#2d5959", "#457ae6", "#99bbff"]
hex_list64 = ["#ffff00", "#a19700", "#b31515", "#457ae6", "#99bbff"]
hex_list65 = ["#ffff00", "#b31515", "#457ae6"]
hex_list66 = ["#1d4e89ff", "#fbd1a2ff", "#7dcfb6ff", "#1d4e89ff"]
hex_list67 = ["#0051f2", "#6da7fb", "#daffa5", "#d9963a", "#b35a00"]
hex_list68 = ["#a00000", "#df7b7b", "#dae695", "#47b0df", "#008fbf"]
hex_list69 = ["#5813fc", "#1cc2fd", "#7dfd94", "#f5c926", "#ff2b18"]
hex_list70 = ["#5813fc", "#1cc2fd", "#4ce600", "#267300", "#f5c926", "#ff2b18"]
hex_list71 = ["#8100e6", "#b360d1", "#f2cf9e", "#6eb830", "#2b9900"]
hex_list72 = ["#5813fc", "#1cc2fd", "#f5c926", "#ff2b18"]
hex_list73 = ["#8100e6", "#b360d1", "#6eb830", "#2b9900"]
hex_list74 = ["#3a7300", "#53a600", "#f09100", "#d957b9", "#ea311f"]
hex_list75 = ["#ffff00", "#a19700", "#b31515", "#457ae6", "#99bbff", "#3a7300"]
hex_list76 = ["#5813fc", "#1cc2fd", "#267300", "#7dfd94", "#f5c926", "#ff2b18"]
hex_list77 = ["9b5de5", "f15bb5", "fee440", "00bbf9", "#e60000"]
hex_list78 = ["9b5de5", "f15bb5", "#e60000", "F9C74F", "00bbf9", "184E77"]
hex_list79 = ["8ecae6", "219ebc", "#e60000", "ffb703", "fb8500", "#54D852"]
hex_list80 = ["#e65154", "#26b6ff", "#67e6d1", "#cd76d6", "#ffca8c", "#fff2b3", "#ff8cd9", "#d99d5b", "#c8f2a9", "#d4b8ff"]
hex_list81 = ["#d92b30", "#0095ba", "#3cccb4", "#ab52b3", "#ffb259", "#ffdf3c", "#eb82eb", "#c27c30", "#a0d17d", "#f260a1"]
hex_list82 = ["#e60049", "#0bb4ff", "#50e991", "#e6d800", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#00bfa0", "#f0cccc"]
hex_list83 = ["#b30000", "#7c1158", "#4421af", "#1a53ff", "#0d88e6", "#00b7c7", "#5ad45a", "#8be04e", "#c5d96d", "#ebdc78"]
hex_list84 = ["04e762", "f5b700", "dc0073", "008bf8", "89fc00"]
hex_list85 = ["ff312e", "119da4", "f5b700", "dc0073", "008bf8", "89fc00", "03071e"]
hex_list86 = ["0d00a4", "f5b700", "dc0073", "008bf8", "89fc00", "011627"]
hex_list87 = ["011627", "f5b700", "dc0073", "0d00a4", "89fc00"]
hex_list88 = ["70e4ef", "#1cc2fd", "#267300", "#7dfd94", "#f5c926", "#ff2b18"]
hex_list89 = ["080708", "3772ff", "df2935", "fdca40", "95c623", "7400b8"]
hex_list90 = ["38b000", "#df7b7b", "#dae695", "#47b0df", "#008fbf", "#e60000"]
hex_list91 = ["38b000", "#df7b7b", "#dae695", "#47b0df", "#e60000"]
hex_list92 = ["9b5de5", "f15bb5", "fee440", "00f5d4", "ff1e35"]
hex_list93 = ["fbaf00","ffd639","ffa3af","007cbe","00af54"]
hex_list94 = ["ffbf00","2274a5","32936f","fc5130","a23b72"]
hex_list95 = ["540d6e","ee4266","ffd23f","3bceac","0ead69"]
hex_list96 = ["5bc0eb","fde74c","9bc53d","e55934","fa7921"]
hex_list97 = ["ffbe0b","fb5607","ff006e","8338ec","3a86ff"]
hex_list98 = ["ffbe0b","fb5607","ff006e","2d00f7","03045e"]
hex_list99 = ["c45f2a","6699cc","fff275","ff8c42","ff3c38"]
hex_list100 = ["ef476f","ffd166","06d6a0","118ab2","073b4c"]
hex_list101 = ["9b5de5","f15bb5","fee440","00bbf9","00f5d4"]
hex_list102 = ["011627","2d00f7","2ec4b6","e71d36","ff9f1c"]
hex_list103 = ["ffbe0b","70e000","ff006e","8338ec","3a86ff"]
hex_list104 = ["ffbe0b","70e000","ff006e","2d00f7","03045e"]
hex_list105 = ["ffbe0b","70e000","ff006e","540d6e","3a86ff"]
hex_list106 = ["10451d","70e000","ff006e","540d6e","3a86ff"]
hex_list107 = ["1a7431","70e000","ff006e","540d6e","3a86ff"]
hex_list108 = ["2c699a", "2c699a", "83e377", "83e377", "83e377",
               "efea5a", "efea5a", "efea5a", "f29e4c", "f29e4c", "f29e4c"]
hex_list109 = ["00ffff", "00ffff", "0000ff", "0000ff", "ff1e35"]
hex_list110 = ["9c4f96", "ff6355", "fba949", "fae442", "8bd448", "2aa8f2"]
hex_list111 = ["ff0017", "ff8900", "ffb200", "ffff00", "94ff00"]
hex_list112 = ["dc8665", "138086", "534666", "cd7672", "eeb462"]
hex_list113 = ["e8a49c", "3c4cad", "240e8b", "240e8b", "f04393", "f04393", "f9c449"]
hex_list114 = ["35bbca", "f8d90f", "f8d90f", "d3dd18", "fe7a15"]
hex_list115 = ["fe7a15", "fe7a15", "f8d90f", "f8d90f", "f8d90f", "f8d90f", "35bbca"]
hex_list116 = ["8ecae6","8ecae6","219ebc","219ebc","023047","023047","023047","fb8500"]
hex_list117 = ["03071e","370617","6a040f","9d0208","d00000","dc2f02","e85d04","f48c06","faa307","ffba08"]
hex_list118 = ["001219","005f73","0a9396","94d2bd","e9d8a6","ee9b00","ca6702","bb3e03","ae2012","9b2226"]
hex_list119 = ["f72585","7209b7","3a0ca3","4361ee","4361ee","4361ee","4cc9f0"]
hex_list120 = ["54478c","2c699a","048ba8","0db39e","16db93","83e377","b9e769","efea5a","f1c453","f29e4c"]
hex_list121 = ["ff499e","d264b6","a480cf","779be7","49b6ff"]
hex_list122 = ["d61c4e", "d61c4e", "fedb39", "fedb39", "1cd6ce", "1cd6ce", "293462"]
hex_list123 = ["89cffd", "fbdf07", "ff7f3f", "f94892"]
hex_list124 = ["ff8b3d", "ff8b3d", "ffda00", "ffda00", "1cd6ce", "1cd6ce", "2a0944"]
hex_list125 = ["6411ad", "6411ad", "6411ad",
               "e3242b", "e3242b", "e3242b",
               "f15bb5", "f15bb5", "f15bb5",
               "fee440", "fee440", "fee440",
               "00bbf9", "00bbf9", "00bbf9"]


