function task_contrasts = get_HCP_task_contrasts

tasks = {'SOCIAL', 'MOTOR', 'GAMBLING', 'WM', 'LANGUAGE', 'EMOTION', 'RELATIONAL'};


%% SOCIAL

task_contrasts.SOCIAL = {'social_random';
                         'social_tom';
                         'social_tom_random'};

%% MOTOR

task_contrasts.MOTOR = {'motor_cue';
                        'motor_lf';
                        'motor_lh';
                        'motor_rf';
                        'motor_rh';
                        'motor_t';
                        'motor_avg';
                        'motor_cue_avg';
                        'motor_lf_avg';
                        'motor_lh_avg';
                        'motor_rf_avg';
                        'motor_rh_avg';
                        'motor_t_avg'};
                    
%% GAMBLING

task_contrasts.GAMBLING = {'gambling_punish';
                           'gambling_reward';
                           'gambling_punish_reward'};
                       
%% WM

task_contrasts.WM = {'wm_2bk_body';
                     'wm_2bk_face';
                     'wm_2bk_place';
                     'wm_2bk_tool';
                     'wm_0bk_body';
                     'wm_0bk_face';
                     'wm_0bk_place';
                     'wm_0bk_tool';
                     'wm_2bk';
                     'wm_0bk';
                     'wm_2bk_0bk';
                     'wm_body';
                     'wm_face';
                     'wm_place';
                     'wm_tool';
                     'wm_body_avg';
                     'wm_face_avg';
                     'wm_place_avg';
                     'wm_tool_avg'};
                       
%% LANGUAGE

task_contrasts.LANGUAGE = {'language_math';
                           'language_story';
                           'language_math_story'};

%% EMOTION

task_contrasts.EMOTION = {'emotion_faces';
                          'emotion_shapes';
                          'emotion_faces_shapes'};
                      
%% RELATIONAL

task_contrasts.RELATIONAL = {'relational_match';
                             'relational_rel';
                             'relational_match_rel'};
                      