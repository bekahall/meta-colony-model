def return_frequency_list(frequency, deme_no):
    if isinstance(frequency, list):
        try:
            for f in frequency:
                if (f < 0) | (f > 1):
                    raise Exception("all freqency values must be between 0 and 1")
        except ValueError:
            print("freqency argument must be a number or a list of numbers")
        return frequency
    else:
        if not isinstance(deme_no, int):
            raise ValueError("deme_no argument must be an integer")
        else:
            try:
                if (frequency < 0) | (frequency > 1):
                    raise Exception("frequency value must be between 0 and 1")
            except ValueError:
                print("frequency argument must be a number of a list of numbers")
            return [frequency] * deme_no
        
def return_cell_count_list(cell_count, deme_no):
    if isinstance(cell_count, list):
        for c in cell_count:
            if not isinstance(c, int):
                raise ValueError("cell_count argument must a positive integer or list of positive integers")
            elif c < 0:
                raise Exception("cell_count argument must a positive integer or list of positive integers")
        return cell_count
    else:
        if not isinstance(deme_no, int):
            raise ValueError("deme_no argument must be an integer")
        else:
            if not isinstance(cell_count, int):
                raise ValueError("cell_count argument must a positive integer or list of positive integers")
            elif cell_count < 0:
                raise Exception("cell_count argument must a positive integer or list of positive integers")
            else:
                return [cell_count] * deme_no
            
def get_deme_no(frequencies, cell_counts):
    freq_length = len(frequencies)
    count_length = len(cell_counts)

    if freq_length != count_length:
        raise Exception("frequency and cell_count must have the same length")
    else:
        return freq_length
    
def check_size(size):
    if not isinstance(size, list):
        raise ValueError("size argument must be a list")
    elif len(size) != 2:
        raise Exception("size argument must have length 2")
    else:
        for s in size:
            if not isinstance(s, int):
                raise ValueError("size argument must contain positive integers")
            elif s < 0:
                raise Exception("size argument must contain positive integers")
        return size