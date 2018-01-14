function [trainedClassifier, validationAccuracy] = trainClassifier(trainingData)
% [trainedClassifier, validationAccuracy] = trainClassifier(trainingData)
% returns a trained classifier and its accuracy. This code recreates the
% classification model trained in Classification Learner app. Use the
% generated code to automate training the same model with new data, or to
% learn how to programmatically train models.
%
%  Input:
%      trainingData: a matrix with the same number of columns and data type
%       as imported into the app.
%
%  Output:
%      trainedClassifier: a struct containing the trained classifier. The
%       struct contains various fields with information about the trained
%       classifier.
%
%      trainedClassifier.predictFcn: a function to make predictions on new
%       data.
%
%      validationAccuracy: a double containing the accuracy in percent. In
%       the app, the History list displays this overall accuracy score for
%       each model.
%
% Use the code to train the model with new data. To retrain your
% classifier, call the function from the command line with your original
% data or new data as the input argument trainingData.
%
% For example, to retrain a classifier trained with the original data set
% T, enter:
%   [trainedClassifier, validationAccuracy] = trainClassifier(T)
%
% To make predictions with the returned 'trainedClassifier' on new data T2,
% use
%   yfit = trainedClassifier.predictFcn(T2)
%
% T2 must be a matrix containing only the predictor columns used for
% training. For details, enter:
%   trainedClassifier.HowToPredict

% Auto-generated by MATLAB on 13-Jan-2018 16:30:37


% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118', 'column_119', 'column_120', 'column_121', 'column_122', 'column_123', 'column_124', 'column_125', 'column_126', 'column_127', 'column_128', 'column_129', 'column_130', 'column_131', 'column_132', 'column_133', 'column_134', 'column_135', 'column_136', 'column_137', 'column_138', 'column_139', 'column_140', 'column_141', 'column_142', 'column_143', 'column_144', 'column_145', 'column_146', 'column_147', 'column_148', 'column_149', 'column_150', 'column_151', 'column_152', 'column_153', 'column_154', 'column_155', 'column_156', 'column_157', 'column_158', 'column_159', 'column_160', 'column_161', 'column_162', 'column_163', 'column_164', 'column_165', 'column_166', 'column_167', 'column_168', 'column_169', 'column_170', 'column_171', 'column_172', 'column_173', 'column_174', 'column_175', 'column_176', 'column_177', 'column_178', 'column_179', 'column_180', 'column_181', 'column_182', 'column_183', 'column_184', 'column_185', 'column_186', 'column_187', 'column_188', 'column_189', 'column_190', 'column_191', 'column_192', 'column_193', 'column_194', 'column_195', 'column_196', 'column_197', 'column_198', 'column_199', 'column_200', 'column_201', 'column_202', 'column_203', 'column_204', 'column_205', 'column_206', 'column_207', 'column_208', 'column_209', 'column_210', 'column_211', 'column_212', 'column_213', 'column_214', 'column_215', 'column_216', 'column_217', 'column_218', 'column_219', 'column_220', 'column_221', 'column_222', 'column_223', 'column_224', 'column_225', 'column_226', 'column_227', 'column_228', 'column_229', 'column_230', 'column_231', 'column_232', 'column_233', 'column_234', 'column_235', 'column_236', 'column_237', 'column_238', 'column_239', 'column_240', 'column_241', 'column_242', 'column_243', 'column_244', 'column_245', 'column_246', 'column_247', 'column_248', 'column_249', 'column_250', 'column_251', 'column_252', 'column_253', 'column_254', 'column_255', 'column_256', 'column_257', 'column_258', 'column_259', 'column_260', 'column_261', 'column_262', 'column_263', 'column_264', 'column_265', 'column_266', 'column_267', 'column_268', 'column_269', 'column_270', 'column_271', 'column_272', 'column_273', 'column_274', 'column_275', 'column_276', 'column_277', 'column_278', 'column_279', 'column_280', 'column_281', 'column_282', 'column_283', 'column_284', 'column_285', 'column_286', 'column_287', 'column_288', 'column_289', 'column_290', 'column_291', 'column_292', 'column_293', 'column_294', 'column_295', 'column_296', 'column_297', 'column_298', 'column_299', 'column_300', 'column_301', 'column_302', 'column_303', 'column_304', 'column_305', 'column_306', 'column_307', 'column_308', 'column_309', 'column_310', 'column_311', 'column_312', 'column_313', 'column_314', 'column_315', 'column_316', 'column_317', 'column_318', 'column_319', 'column_320', 'column_321', 'column_322', 'column_323', 'column_324', 'column_325', 'column_326', 'column_327', 'column_328', 'column_329', 'column_330', 'column_331', 'column_332', 'column_333', 'column_334', 'column_335', 'column_336', 'column_337', 'column_338', 'column_339', 'column_340', 'column_341', 'column_342', 'column_343', 'column_344', 'column_345', 'column_346', 'column_347', 'column_348', 'column_349', 'column_350', 'column_351', 'column_352', 'column_353', 'column_354', 'column_355', 'column_356', 'column_357', 'column_358', 'column_359', 'column_360', 'column_361', 'column_362', 'column_363', 'column_364', 'column_365', 'column_366', 'column_367', 'column_368', 'column_369', 'column_370', 'column_371', 'column_372', 'column_373', 'column_374', 'column_375', 'column_376', 'column_377', 'column_378', 'column_379', 'column_380', 'column_381', 'column_382', 'column_383', 'column_384', 'column_385', 'column_386', 'column_387', 'column_388', 'column_389', 'column_390', 'column_391', 'column_392', 'column_393', 'column_394', 'column_395', 'column_396', 'column_397', 'column_398', 'column_399', 'column_400', 'column_401', 'column_402', 'column_403', 'column_404', 'column_405', 'column_406', 'column_407', 'column_408', 'column_409', 'column_410', 'column_411', 'column_412', 'column_413', 'column_414', 'column_415', 'column_416', 'column_417', 'column_418', 'column_419', 'column_420', 'column_421', 'column_422', 'column_423', 'column_424', 'column_425', 'column_426', 'column_427', 'column_428', 'column_429', 'column_430', 'column_431', 'column_432', 'column_433', 'column_434', 'column_435', 'column_436', 'column_437', 'column_438', 'column_439', 'column_440', 'column_441', 'column_442'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118', 'column_119', 'column_120', 'column_121', 'column_122', 'column_123', 'column_124', 'column_125', 'column_126', 'column_127', 'column_128', 'column_129', 'column_130', 'column_131', 'column_132', 'column_133', 'column_134', 'column_135', 'column_136', 'column_137', 'column_138', 'column_139', 'column_140', 'column_141', 'column_142', 'column_143', 'column_144', 'column_145', 'column_146', 'column_147', 'column_148', 'column_149', 'column_150', 'column_151', 'column_152', 'column_153', 'column_154', 'column_155', 'column_156', 'column_157', 'column_158', 'column_159', 'column_160', 'column_161', 'column_162', 'column_163', 'column_164', 'column_165', 'column_166', 'column_167', 'column_168', 'column_169', 'column_170', 'column_171', 'column_172', 'column_173', 'column_174', 'column_175', 'column_176', 'column_177', 'column_178', 'column_179', 'column_180', 'column_181', 'column_182', 'column_183', 'column_184', 'column_185', 'column_186', 'column_187', 'column_188', 'column_189', 'column_190', 'column_191', 'column_192', 'column_193', 'column_194', 'column_195', 'column_196', 'column_197', 'column_198', 'column_199', 'column_200', 'column_201', 'column_202', 'column_203', 'column_204', 'column_205', 'column_206', 'column_207', 'column_208', 'column_209', 'column_210', 'column_211', 'column_212', 'column_213', 'column_214', 'column_215', 'column_216', 'column_217', 'column_218', 'column_219', 'column_220', 'column_221', 'column_222', 'column_223', 'column_224', 'column_225', 'column_226', 'column_227', 'column_228', 'column_229', 'column_230', 'column_231', 'column_232', 'column_233', 'column_234', 'column_235', 'column_236', 'column_237', 'column_238', 'column_239', 'column_240', 'column_241', 'column_242', 'column_243', 'column_244', 'column_245', 'column_246', 'column_247', 'column_248', 'column_249', 'column_250', 'column_251', 'column_252', 'column_253', 'column_254', 'column_255', 'column_256', 'column_257', 'column_258', 'column_259', 'column_260', 'column_261', 'column_262', 'column_263', 'column_264', 'column_265', 'column_266', 'column_267', 'column_268', 'column_269', 'column_270', 'column_271', 'column_272', 'column_273', 'column_274', 'column_275', 'column_276', 'column_277', 'column_278', 'column_279', 'column_280', 'column_281', 'column_282', 'column_283', 'column_284', 'column_285', 'column_286', 'column_287', 'column_288', 'column_289', 'column_290', 'column_291', 'column_292', 'column_293', 'column_294', 'column_295', 'column_296', 'column_297', 'column_298', 'column_299', 'column_300', 'column_301', 'column_302', 'column_303', 'column_304', 'column_305', 'column_306', 'column_307', 'column_308', 'column_309', 'column_310', 'column_311', 'column_312', 'column_313', 'column_314', 'column_315', 'column_316', 'column_317', 'column_318', 'column_319', 'column_320', 'column_321', 'column_322', 'column_323', 'column_324', 'column_325', 'column_326', 'column_327', 'column_328', 'column_329', 'column_330', 'column_331', 'column_332', 'column_333', 'column_334', 'column_335', 'column_336', 'column_337', 'column_338', 'column_339', 'column_340', 'column_341', 'column_342', 'column_343', 'column_344', 'column_345', 'column_346', 'column_347', 'column_348', 'column_349', 'column_350', 'column_351', 'column_352', 'column_353', 'column_354', 'column_355', 'column_356', 'column_357', 'column_358', 'column_359', 'column_360', 'column_361', 'column_362', 'column_363', 'column_364', 'column_365', 'column_366', 'column_367', 'column_368', 'column_369', 'column_370', 'column_371', 'column_372', 'column_373', 'column_374', 'column_375', 'column_376', 'column_377', 'column_378', 'column_379', 'column_380', 'column_381', 'column_382', 'column_383', 'column_384', 'column_385', 'column_386', 'column_387', 'column_388', 'column_389', 'column_390', 'column_391', 'column_392', 'column_393', 'column_394', 'column_395', 'column_396', 'column_397', 'column_398', 'column_399', 'column_400', 'column_401', 'column_402', 'column_403', 'column_404', 'column_405', 'column_406', 'column_407', 'column_408', 'column_409', 'column_410', 'column_411', 'column_412', 'column_413', 'column_414', 'column_415', 'column_416', 'column_417', 'column_418', 'column_419', 'column_420', 'column_421', 'column_422', 'column_423', 'column_424', 'column_425', 'column_426', 'column_427', 'column_428', 'column_429', 'column_430', 'column_431', 'column_432', 'column_433', 'column_434', 'column_435', 'column_436', 'column_437', 'column_438', 'column_439', 'column_440', 'column_441'};
predictors = inputTable(:, predictorNames);
response = inputTable.column_442;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
template = templateSVM(...
    'KernelFunction', 'linear', ...
    'PolynomialOrder', [], ...
    'KernelScale', 'auto', ...
    'BoxConstraint', 1, ...
    'Standardize', true);
classificationSVM = fitcecoc(...
    predictors, ...
    response, ...
    'Learners', template, ...
    'Coding', 'onevsone', ...
    'ClassNames', [6; 7; 8; 9; 12; 14; 17; 18; 20; 21; 22; 23; 24; 25; 26; 27; 28; 29; 30; 31; 32; 33; 34; 35; 37; 38; 39; 40; 41; 43; 44; 45; 46; 47; 48; 49; 50; 52; 53; 55; 56; 57; 58; 59; 60; 61; 62; 63; 64; 66; 67; 68; 69; 70; 71; 72; 73; 74; 76; 77; 78; 79; 80; 83; 84; 85; 86; 87; 88; 89; 90; 91; 92]);

% Create the result struct with predict function
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
svmPredictFcn = @(x) predict(classificationSVM, x);
trainedClassifier.predictFcn = @(x) svmPredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.ClassificationSVM = classificationSVM;
trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2017b.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 441 columns because this model was trained using 441 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118', 'column_119', 'column_120', 'column_121', 'column_122', 'column_123', 'column_124', 'column_125', 'column_126', 'column_127', 'column_128', 'column_129', 'column_130', 'column_131', 'column_132', 'column_133', 'column_134', 'column_135', 'column_136', 'column_137', 'column_138', 'column_139', 'column_140', 'column_141', 'column_142', 'column_143', 'column_144', 'column_145', 'column_146', 'column_147', 'column_148', 'column_149', 'column_150', 'column_151', 'column_152', 'column_153', 'column_154', 'column_155', 'column_156', 'column_157', 'column_158', 'column_159', 'column_160', 'column_161', 'column_162', 'column_163', 'column_164', 'column_165', 'column_166', 'column_167', 'column_168', 'column_169', 'column_170', 'column_171', 'column_172', 'column_173', 'column_174', 'column_175', 'column_176', 'column_177', 'column_178', 'column_179', 'column_180', 'column_181', 'column_182', 'column_183', 'column_184', 'column_185', 'column_186', 'column_187', 'column_188', 'column_189', 'column_190', 'column_191', 'column_192', 'column_193', 'column_194', 'column_195', 'column_196', 'column_197', 'column_198', 'column_199', 'column_200', 'column_201', 'column_202', 'column_203', 'column_204', 'column_205', 'column_206', 'column_207', 'column_208', 'column_209', 'column_210', 'column_211', 'column_212', 'column_213', 'column_214', 'column_215', 'column_216', 'column_217', 'column_218', 'column_219', 'column_220', 'column_221', 'column_222', 'column_223', 'column_224', 'column_225', 'column_226', 'column_227', 'column_228', 'column_229', 'column_230', 'column_231', 'column_232', 'column_233', 'column_234', 'column_235', 'column_236', 'column_237', 'column_238', 'column_239', 'column_240', 'column_241', 'column_242', 'column_243', 'column_244', 'column_245', 'column_246', 'column_247', 'column_248', 'column_249', 'column_250', 'column_251', 'column_252', 'column_253', 'column_254', 'column_255', 'column_256', 'column_257', 'column_258', 'column_259', 'column_260', 'column_261', 'column_262', 'column_263', 'column_264', 'column_265', 'column_266', 'column_267', 'column_268', 'column_269', 'column_270', 'column_271', 'column_272', 'column_273', 'column_274', 'column_275', 'column_276', 'column_277', 'column_278', 'column_279', 'column_280', 'column_281', 'column_282', 'column_283', 'column_284', 'column_285', 'column_286', 'column_287', 'column_288', 'column_289', 'column_290', 'column_291', 'column_292', 'column_293', 'column_294', 'column_295', 'column_296', 'column_297', 'column_298', 'column_299', 'column_300', 'column_301', 'column_302', 'column_303', 'column_304', 'column_305', 'column_306', 'column_307', 'column_308', 'column_309', 'column_310', 'column_311', 'column_312', 'column_313', 'column_314', 'column_315', 'column_316', 'column_317', 'column_318', 'column_319', 'column_320', 'column_321', 'column_322', 'column_323', 'column_324', 'column_325', 'column_326', 'column_327', 'column_328', 'column_329', 'column_330', 'column_331', 'column_332', 'column_333', 'column_334', 'column_335', 'column_336', 'column_337', 'column_338', 'column_339', 'column_340', 'column_341', 'column_342', 'column_343', 'column_344', 'column_345', 'column_346', 'column_347', 'column_348', 'column_349', 'column_350', 'column_351', 'column_352', 'column_353', 'column_354', 'column_355', 'column_356', 'column_357', 'column_358', 'column_359', 'column_360', 'column_361', 'column_362', 'column_363', 'column_364', 'column_365', 'column_366', 'column_367', 'column_368', 'column_369', 'column_370', 'column_371', 'column_372', 'column_373', 'column_374', 'column_375', 'column_376', 'column_377', 'column_378', 'column_379', 'column_380', 'column_381', 'column_382', 'column_383', 'column_384', 'column_385', 'column_386', 'column_387', 'column_388', 'column_389', 'column_390', 'column_391', 'column_392', 'column_393', 'column_394', 'column_395', 'column_396', 'column_397', 'column_398', 'column_399', 'column_400', 'column_401', 'column_402', 'column_403', 'column_404', 'column_405', 'column_406', 'column_407', 'column_408', 'column_409', 'column_410', 'column_411', 'column_412', 'column_413', 'column_414', 'column_415', 'column_416', 'column_417', 'column_418', 'column_419', 'column_420', 'column_421', 'column_422', 'column_423', 'column_424', 'column_425', 'column_426', 'column_427', 'column_428', 'column_429', 'column_430', 'column_431', 'column_432', 'column_433', 'column_434', 'column_435', 'column_436', 'column_437', 'column_438', 'column_439', 'column_440', 'column_441', 'column_442'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118', 'column_119', 'column_120', 'column_121', 'column_122', 'column_123', 'column_124', 'column_125', 'column_126', 'column_127', 'column_128', 'column_129', 'column_130', 'column_131', 'column_132', 'column_133', 'column_134', 'column_135', 'column_136', 'column_137', 'column_138', 'column_139', 'column_140', 'column_141', 'column_142', 'column_143', 'column_144', 'column_145', 'column_146', 'column_147', 'column_148', 'column_149', 'column_150', 'column_151', 'column_152', 'column_153', 'column_154', 'column_155', 'column_156', 'column_157', 'column_158', 'column_159', 'column_160', 'column_161', 'column_162', 'column_163', 'column_164', 'column_165', 'column_166', 'column_167', 'column_168', 'column_169', 'column_170', 'column_171', 'column_172', 'column_173', 'column_174', 'column_175', 'column_176', 'column_177', 'column_178', 'column_179', 'column_180', 'column_181', 'column_182', 'column_183', 'column_184', 'column_185', 'column_186', 'column_187', 'column_188', 'column_189', 'column_190', 'column_191', 'column_192', 'column_193', 'column_194', 'column_195', 'column_196', 'column_197', 'column_198', 'column_199', 'column_200', 'column_201', 'column_202', 'column_203', 'column_204', 'column_205', 'column_206', 'column_207', 'column_208', 'column_209', 'column_210', 'column_211', 'column_212', 'column_213', 'column_214', 'column_215', 'column_216', 'column_217', 'column_218', 'column_219', 'column_220', 'column_221', 'column_222', 'column_223', 'column_224', 'column_225', 'column_226', 'column_227', 'column_228', 'column_229', 'column_230', 'column_231', 'column_232', 'column_233', 'column_234', 'column_235', 'column_236', 'column_237', 'column_238', 'column_239', 'column_240', 'column_241', 'column_242', 'column_243', 'column_244', 'column_245', 'column_246', 'column_247', 'column_248', 'column_249', 'column_250', 'column_251', 'column_252', 'column_253', 'column_254', 'column_255', 'column_256', 'column_257', 'column_258', 'column_259', 'column_260', 'column_261', 'column_262', 'column_263', 'column_264', 'column_265', 'column_266', 'column_267', 'column_268', 'column_269', 'column_270', 'column_271', 'column_272', 'column_273', 'column_274', 'column_275', 'column_276', 'column_277', 'column_278', 'column_279', 'column_280', 'column_281', 'column_282', 'column_283', 'column_284', 'column_285', 'column_286', 'column_287', 'column_288', 'column_289', 'column_290', 'column_291', 'column_292', 'column_293', 'column_294', 'column_295', 'column_296', 'column_297', 'column_298', 'column_299', 'column_300', 'column_301', 'column_302', 'column_303', 'column_304', 'column_305', 'column_306', 'column_307', 'column_308', 'column_309', 'column_310', 'column_311', 'column_312', 'column_313', 'column_314', 'column_315', 'column_316', 'column_317', 'column_318', 'column_319', 'column_320', 'column_321', 'column_322', 'column_323', 'column_324', 'column_325', 'column_326', 'column_327', 'column_328', 'column_329', 'column_330', 'column_331', 'column_332', 'column_333', 'column_334', 'column_335', 'column_336', 'column_337', 'column_338', 'column_339', 'column_340', 'column_341', 'column_342', 'column_343', 'column_344', 'column_345', 'column_346', 'column_347', 'column_348', 'column_349', 'column_350', 'column_351', 'column_352', 'column_353', 'column_354', 'column_355', 'column_356', 'column_357', 'column_358', 'column_359', 'column_360', 'column_361', 'column_362', 'column_363', 'column_364', 'column_365', 'column_366', 'column_367', 'column_368', 'column_369', 'column_370', 'column_371', 'column_372', 'column_373', 'column_374', 'column_375', 'column_376', 'column_377', 'column_378', 'column_379', 'column_380', 'column_381', 'column_382', 'column_383', 'column_384', 'column_385', 'column_386', 'column_387', 'column_388', 'column_389', 'column_390', 'column_391', 'column_392', 'column_393', 'column_394', 'column_395', 'column_396', 'column_397', 'column_398', 'column_399', 'column_400', 'column_401', 'column_402', 'column_403', 'column_404', 'column_405', 'column_406', 'column_407', 'column_408', 'column_409', 'column_410', 'column_411', 'column_412', 'column_413', 'column_414', 'column_415', 'column_416', 'column_417', 'column_418', 'column_419', 'column_420', 'column_421', 'column_422', 'column_423', 'column_424', 'column_425', 'column_426', 'column_427', 'column_428', 'column_429', 'column_430', 'column_431', 'column_432', 'column_433', 'column_434', 'column_435', 'column_436', 'column_437', 'column_438', 'column_439', 'column_440', 'column_441'};
predictors = inputTable(:, predictorNames);
response = inputTable.column_442;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% Perform cross-validation
partitionedModel = crossval(trainedClassifier.ClassificationSVM, 'KFold', 5);

% Compute validation predictions
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');
