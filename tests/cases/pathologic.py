import mpwt.mpwt_workflow
from multiprocessing import freeze_support

# print(mpwt.multiprocess_pwt(input_folder='pathologic_template',
# patho_inference=True))
if __name__ == "__main__":
    freeze_support()
    mpwt.run_pwt("pathologic_template", False, False, False, False, True, True)
