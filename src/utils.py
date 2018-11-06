from enum import Enum


class InputImageType(Enum):
    Tiff = 0
    OmeTiff = 1
    ImageJTiff = 2
    MultipleTiff = 3
    Nifti = 4


class ImageComponentType(Enum):
    SCALAR = 0
    VECTOR = 1
    SEGMENTATION = 2


class PathUtil(object):

    _group_format_str = 't{:0>5}/s{:0>2}/{:}'
    _cells_format_str = '{:}/cells'
    _res_format_str = 's{:0>2}/resolutions'
    _sub_format_str = 's{:0>2}/subdivisions'
    _affine_format_str = 's{:0>2}/affine'
    _setup_format_str = 's{:0>2}'
    _lsfm_id_format_str = 'LSFM/{:}'
    _lsfm_image_format_str = 'LSFM/{:}/IMAGE'
    _lsfm_image_level_str = 'LSFM/{:}/IMAGE/{:}/cells'
    _lsfm_metadata_format_str = 'LSFM/{:}/METADATA'
    _lsfm_affine_format_str = 'LSFM/{:}/AFFINES/{:}'
    _lsfm_affines_group_format_str = 'LSFM/{:}/AFFINES'
    _lsfm_affine_itk_params_format_str = 'LSFM/{:}/AFFINES/ITK_PARAMS/{:}'
    _lsfm_affine_itk_fixed_params_format_str = 'LSFM/{:}/AFFINES/FIXED_ITK_PARAMS/{:}'
    _lsfm_deformation_field_format_str = 'LSFM/{:}/DEFORMATION_FIELDS/{:}'
    _lsfm_deformation_fields_str = 'LSFM/{:}/DEFORMATION_FIELDS'
    _lsfm_deformation_field_data_format_str = 'LSFM/{:}/DEFORMATION_FIELDS/{:}/cells'
    _log_format_str = 'log/{:%Y-%m-%d_%H:%M:%S}_{:}'
    _logs_str = 'log/'
    _extended_metadata_str = '/extended_metadata'

    @classmethod
    def get_extended_metadata_path(cls):
        return cls._extended_metadata_str

    @classmethod
    def get_logs_path(cls):
        return cls._logs_str

    @classmethod
    def get_log_path(cls, timestamp, cmd_name):
        return cls._log_format_str.format(timestamp,
                                          cmd_name)

    @classmethod
    def get_lsfm_id_path(cls, _id):
        return cls._lsfm_id_format_str.format(_id)

    @classmethod
    def get_setup_path(cls, setup):
        return cls._setup_format_str.format(setup)

    @classmethod
    def get_group_path(cls, timepoint, setup, level):
        return cls._group_format_str.format(timepoint,
                                            setup,
                                            level)

    @classmethod
    def get_cells_path(cls, timepoint, setup, level):
        return cls._cells_format_str.format(cls.get_group_path(timepoint,
                                                               setup,
                                                               level))

    @classmethod
    def get_res_path(cls, setup):
        return cls._res_format_str.format(setup)

    @classmethod
    def get_sub_path(cls, setup):
        return cls._sub_format_str.format(setup)

    @classmethod
    def get_affine_path(cls, setup):
        return cls._affine_format_str.format(setup)

    @classmethod
    def get_lsfm_image_path(cls, case_id):
        return cls._lsfm_image_format_str.format(case_id)

    @classmethod
    def get_lsfm_image_cells_path(cls, case_id, level):
        return cls._lsfm_image_level_str.format(case_id, level)

    @classmethod
    def get_lsfm_affine_path(cls, case_id, affine_name):
        return cls._lsfm_affine_format_str.format(case_id, affine_name)

    @classmethod
    def get_lsfm_affines_group_path(cls, case_id):
        return cls._lsfm_affines_group_format_str.format(case_id)

    @classmethod
    def get_lsfm_affine_itk_path(cls, case_id, affine_name):
        return (cls._lsfm_affine_itk_params_format_str.format(case_id, affine_name),
                cls._lsfm_affine_itk_fixed_params_format_str.format(case_id, affine_name))

    @classmethod
    def get_lsfm_df_path(cls, case_id, df_name):
        return cls._lsfm_deformation_field_format_str.format(case_id, df_name)

    @classmethod
    def get_lsfm_df_data_path(cls, case_id, df_name):
        return cls._lsfm_deformation_field_data_format_str.format(case_id, df_name)



