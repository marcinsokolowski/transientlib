/*********************  ERROR codes **********************/
#define E_no_error		0
#define E_OK			0
#define E_ok			0
#define E_timer			1
#define E_invalid_motor		2
#define E_invalid_pulse		3
#define E_invalid_offset	4
#define E_invalid_position	5
#define E_object_too_low	6
#define E_motor_busy		7
#define E_invalid_accel		8
#define E_invalid_pmax		9
#define E_invalid_pmin		10
#define E_invalid_encoder	11
#define E_cam_not_connected	12
#define E_mot_not_connected	13
#define E_open_cam		14
#define E_open_mot		15
#define E_baud_failed		16
#define E_illegal_sbinn		17/* illegal software binning specified */
#define E_illegal_abinn		18/* illegal software binning specified */
#define E_parameters		19/* bad nuber of parameters */

#define E_camera_busy		21	/* camera is busy - do not disturb it */
#define E_open_macro		22	/* macro file does not exist */
#define E_cannot_open_serial    23
#define E_cannot_continue_macro	24
#define E_missing_arg		25
#define E_illegal_arg		26
#define E_limit_switch		27
#define E_mot_not_responding	28
#define E_cam_not_responding	29
#define E_db_problem		30
#define E_soft_problem		31
#define E_hard_problem		32
#define E_sun_too_high		33
#define E_moon_too_close	34
#define E_mask_obscuration	35
#define E_sync_timeout		36
#define E_sync_motor		37
#define E_aborted		38
#define E_malloc		39
#define E_invalid_port          40
#define E_invalid_value         41


#define E_create_fits  		42
#define E_no_dark_correction	43
#define E_no_flat_correction 	44
#define E_no_dark       	45
#define E_no_bias       	46
#define E_no_flat       	47
#define E_dark_done     	48
#define E_flat_done     	49
#define E_psf_not_odd   	50
#define E_diff_size		51
#define E_missing_keyword	52
#define E_no_cat_stars		53
#define E_no_common_stars	54
#define E_cannot_transform	55
#define E_no_stars		56
#define E_imroper_frame		57
#define E_no_photometry		58
#define E_photometry_done	59
#define E_no_astrometry		60
#define E_astrometry_done	61
#define E_open_fits		62
#define E_read_fits  		63	
#define E_saturated		64
#define E_open_file		65
#define E_grb			67		/* GRB event arrived */

#define E_environment		68

/* sky_level output codes */
#define E_sky_error		70
#define E_sky_toolow		71
#define E_sky_toohigh		72
#define E_sky_toolong		73
#define E_sky_tooshort		74

/* database erros */
#define E_open_flog		100
#define E_cannot_add_flog	101
#define E_already_in_flog	102

/* astrometry timeout */
#define E_astrometry_timeout 103

/* to small number of matches - when seems that converged but small number */
#define E_num_matches_too_small 104

/* */
#define E_stop_forced_externally 105

/* timeouts */
#define T_NONE		0x00000000
#define T_MSEND		0x00000001
#define T_MOTASK	0x00000002
#define T_TRY		0x00000004
#define T_BAUD_RATE	0x00000008
#define T_QUIT		0x00000010
#define T_CONFIRM	0x00000020
#define T_CONFIRM2	0x00000040
#define T_NEW_BAUD_RATE	0x00000080
#define T_DOME_OPEN     0x00000100
#define T_DOME_CLOSE    0x00000200

