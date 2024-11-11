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
hex_list126 = ["f94144","f3722c","f8961e","f9c74f","90be6d","43aa8b","577590"]
hex_list127 = ["f72585","b5179e","7209b7","560bad","480ca8","3a0ca3","3f37c9","4361ee","4895ef","4cc9f0"]
hex_list128 = ["ffbe0b","fb5607","ff006e","8338ec","3a86ff"]
hex_list129 = ["f72585","f72585","7209b7","3a0ca3","4361ee","4cc9f0"]
hex_list130 = ["b7094c","a01a58","892b64","723c70","5c4d7d","455e89","2e6f95","1780a1","0091ad"]
hex_list131 = ["54478c","2c699a","048ba8","0db39e","16db93","83e377","b9e769","efea5a","f1c453","f29e4c"]
hex_list132 = ['#4aad5a', '#b5d663', '#ffe710', '#ff9c08', '#ff9c08', '#ff7b10', '#ff7b10']
hex_list133 = ["247ba0","70c1b3","b2dbbf","f3ffbd","ff1654"]
hex_list134 = ["8ecae6","219ebc","126782","023047","ffb703","fd9e02","fb8500"]
hex_list135 = ["042a2b","5eb1bf","cdedf6","ef7b45","d84727"]
hex_list136 = ['#440154', '#3b528b', '#21918c', '#5ec962', '#fde725', '#fde725', '#fde725']
hex_list137 = ["8ecae6","219ebc","126782","023047","ffb703","fd9e02","fb8500"]
hex_list138 = ["ff499e","d264b6","a480cf","779be7","49b6ff"]
hex_list139 = ["247ba0","247ba0","70c1b3","b2dbbf","f3ffbd","ff1654"]
hex_list140 = ['#0a0079', '#0d81f8', '#44caff', '#8aecae', '#8aecae', '#dff58d', '#dff58d','#dff58d','#ffa045',
               '#ffa045', '#ffa045','#ff5a5a', '#ff5a5a', '#ff5a5a']
hex_list141 = ["f79256","fbd1a2","7dcfb6","00b2ca","1d4e89"]
hex_list142 = ['#9e0142','#fee08b','#e6f598','#abdda4','#3288bd','#3288bd','#3288bd']
hex_list143 = ['#440154', '#472878',
               '#1e9e89', '#6cce59', '#6cce59','#b4de2c', '#b4de2c','#fde725','#fde725','#fde725']
hex_list144 = ['#0b0405', '#3e356b', '#357ba3', '#4bc2ad', '#c9decf']
hex_list145 = ['#e41a1c', '#e41a1c', '#e41a1c', '#377eb8', '#377eb8', '#984ea3', '#984ea3', '#ff7f00', '#ff7f00',
               '#ebeb2f']
hex_list146 = ["083d77","ebebd3","f4d35e","ee964b","f95738"]
hex_list147 = ["ff595e","ffca3a","8ac926","1982c4","6a4c93"]
hex_list148 = ['#0000ff','#0000ff','#0000ff','#0000ff',
               '#0059ff',
               '#ffff00','#ffff00','#ffff00',
               '#ff8000',
               '#e30000']
hex_list149 = ['#0a0079', '#0d81f8', '#dff58d', '#dff58d','#dff58d','#ffa045',
               '#ffa045', '#ffa045','#ff5a5a', '#ff5a5a', '#ff5a5a']
hex_list150 = ['#000000', '#780089', '#7a009b', '#0000b2', '#0019dd',
               '#0080dd', '#009ecd', '#00aa9d', '#00a34f', '#00a900',
               '#00cd00', '#00f100', '#76ff00', '#def300', '#fbd500',
               '#ffa400', '#ff1800', '#e10000', '#cd0000', '#cccccc']
hex_list151 = ['#920000', '#f2ff00', '#05f80d', '#3ab594', '#5973ff',
               '#252fff', '#1400ff', '#5800ff', '#9b00ff', '#de00ff',
               '#ff00f2', '#ff00d7', '#ff00bc', '#ff00a1', '#ff0086',
               '#ff006b', '#ff0051', '#ff0036', '#ff001b', '#ff0000']
hex_list152 = ['#640065', '#8200bc', '#5c00ff', '#0021ff', '#0092ff',
               '#00eeff', '#00ff37', '#56ff00', '#9fff00', '#e1ff00',
               '#ffdb00', '#ff9400', '#ff4200', '#ff0000', '#ff0000',
               '#ff0000', '#e00000', '#b80000', '#8e0000', '#610000']
hex_list153 = ['#400040', '#4000c0', '#0040ff', '#0080ff', '#00a0ff',
                '#40c0ff', '#40e0ff', '#40ffff', '#40ffc0', '#40ff40',
                '#80ff40', '#c0ff40', '#ffff40', '#ffe040', '#ffa040',
                '#ff6040', '#ff2040', '#ff60c0', '#ffa0ff', '#ffe0e1']
hex_list154 = ['#141414', '#41301f', '#6e4c2b', '#9b6836', '#c27f3c',
                '#d37f2b', '#e47f1b', '#f57f0a', '#ff8c00', '#ffae00',
                '#ffd000', '#fff200', '#d7ff00', '#94ff00', '#51ff00',
                '#0dff00', '#00f228', '#00e15a', '#00d08d', '#00bfbf']
hex_list155 = ['#3d31ff', '#ff80fc', '#fac1b9', '#f4f288', '#68fa8c',
                '#14c391', '#1d6d94', '#655db8', '#d970f2', '#ed68e5',
                '#ef5ccb', '#f14fb0', '#f44296', '#f6377f', '#f72e6c',
                '#f92559', '#fa1b46', '#fc1234', '#fd0921', '#ff000e']
hex_list156 = ['#35b779', '#654a06', '#f7b719', '#fc7502', '#f44905', #00ff00
                '#e10b11', '#bc0c29', '#970d41', '#790e55', '#620f64',
                '#4b1072', '#341081', '#1d1190', '#1c1191', '#1c118f',
                '#1b108c', '#1a1089', '#1a1084', '#180e7d', '#0000ff']
hex_list157 = ['#28dce2', '#52a9dd', '#7b75d7', '#8f5dd5', '#8f5dd5',
                '#915fd4', '#ab86bf', '#c6adab', '#e1d497', '#fcfa82',
                '#facd63', '#f59541', '#e17b44', '#c96950', '#b1565d',
                '#994469', '#893871', '#792c79', '#681f81', '#581389']
hex_list158 = ['#0000ff', '#0028d7', '#0051ae', '#007986', '#00a15e',
                '#00c936', '#00f20d', '#1bff00', '#43ff00', '#6bff00',
                '#94ff00', '#bcff00', '#e4ff00', '#fff200', '#ffc900',
                '#ffa100', '#ff7900', '#ff5100', '#ff2800', '#ff0000']
hex_list159 = ['#ffff00', '#ffc700', '#ff8e00', '#ff5500', '#ff1c00',
                '#e30000', '#aa0000', '#720000', '#390000', '#000000']
hex_list160 = ['#0000ff', '#1700c1', '#2d0084', '#440046', '#5a0008',
                '#9f0000', '#d83c16', '#f77a2c', '#f77a2c', '#f3d700'] # #f9e317 #fff783
hex_list161 = ['#0000ff', '#3900c6', '#71008e', '#aa0055', '#e3001c',
                '#ff1c00', '#ff5500', '#ff8e00', '#ffc600', '#e9e900']

hex_list162 =  ['#00ff31', '#654a06', '#f7b719', '#fc7502', '#f44905',
                '#e10b11', '#bc0c29', '#970d41', '#790e55', '#620f64',
                '#4b1072', '#341081', '#1d1190', '#1c1191', '#1c118f',
                '#1b108c', '#1a1089', '#1a1084', '#180e7d', '#0000ff']

hex_list163 = ['#5e4fa2', '#4b67ad', '#397fb9', '#4097b7', '#56b0ad',
               '#6dc5a5', '#8ad0a4', '#a7dca4', '#c1e6a0', '#daf09b',
               '#fed582', '#fdc070', '#fdab5f', '#f98f53', '#f57446',
               '#ea5e47', '#dd4a4c', '#cc344d', '#b51b47', '#9e0142']

hex_list164 = ['#ff006f', '#bf3d93', '#8079b7', '#40b6db', '#00f2ff']

hex_list165 = ['#f0f921', '#f7e225', '#fccb26', '#feb72d', '#fba338',
               '#f69044', '#f07e4f', '#e76d5b', '#dc5e66', '#d14e72',
               '#c53f7e', '#b7308b', '#a72197', '#9612a1', '#8305a7',
               '#6f00a8', '#5901a5', '#44049e', '#2c0594', '#0d0887']

hex_list166 = ['#00ff00', '#39ff00', '#71ff00', '#aaff00', '#e3ff00',
               '#ffe300', '#ffaa00' , '#ff7100', '#ff3900', '#ff0000']

hex_list167 = ['#f6f620', '#fdd523', '#ff00ff', '#8080ff', '#00ffff']

hex_list168 = ['#0000ff', '#0000ff', '#ffff00', '#ffff00', '#ff0000']

hex_list169 = ['#ff0000', '#ffbf00', '#80ff7f', '#00bfff', '#0000ff']

hex_list170 = ['#940a02', '#940a02', '#ff00ff', '#ff00ff', '#7f80ff', '#00ffff']

hex_list171 = ['#940a02', '#00ffff']

hex_list172 = ["f94144","f3722c","f8961e","f9c74f","90be6d","43aa8b","577590"]

hex_list173 = ["ff595e","ffca3a","8ac926","1982c4","6a4c93"]

hex_list174 = ["001219","005f73","0a9396","94d2bd","e9d8a6","ee9b00","ca6702","bb3e03","ae2012","9b2226"]

hex_list175 = ["390099","9e0059","ff0054","ff5400","ffbd00"]

hex_list176 = ["f72585","b5179e","7209b7","560bad","480ca8","3a0ca3","3f37c9","4361ee","4895ef","4cc9f0"]

hex_list177 = ["f72585","b5179e","7209b7","560bad","480ca8","3a0ca3","3f37c9","4361ee","4895ef","4cc9f0"]

hex_list178 = ["247ba0","70c1b3","b2dbbf","f3ffbd","ff1654"]

hex_list179 = ["ffbe0b","fb5607","ff006e","8338ec","3a86ff"]

hex_list180 = ["023047","023047","023047","023047","219ebc","219ebc","ffb703"]

hex_list181 = ["ffbe0b","ffbe0b","fb5607","fb5607","ff006e","ff006e","8338ec","4cc9f0"]

hex_list182 = ["ffbe0b","fb5607","ff006e","8338ec","3a0ca3"]




















