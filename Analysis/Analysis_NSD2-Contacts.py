from __future__ import print_function
import mdtraj as mdt
import os
import re
import numpy as np
import collections
from collections import defaultdict
import pandas as pd
import itertools
from contact_map import ContactMap, ContactFrequency, ContactDifference, ResidueContactConcurrence, plot_concurrence
import glob
import pickle
import warnings

warnings.filterwarnings("ignore", message="top= kwargs ignored since this file parser does not support it")

# Input parameters. Needed to find and name the files.

Protein = 'NSD2'        # name of the protein.
Peptide = 'H3K36'       # Name of the complexed peptide
sim_time = '100ns'

number_replicates = 1 # number of replicates to analyze

# all trajectories in the target folder will be loaded, joined and superposed to one big trajectory, which is then analyzed
traj_dict = {}
for i in range(number_replicates):
    folder = '/path-to-trajectory/production_{}_{}_{}_{}.h5'.format(Protein, Peptide, sim_time, i+1)
    print(folder)
    traj_dict[i+1]=mdt.load(folder)

traj_list = []
for key in traj_dict:
    traj_list.append(traj_dict[key])

traj_pre = mdt.join(traj_list,check_topology=True, discard_overlapping_frames=True)
traj = traj_pre.superpose(traj_pre,frame=0,parallel=True)
topology = traj.topology
print('Trajectories successfully loaded, joined and superposed')

# select peptide and protein
protein_peptide1 = topology.select('chainid 0 or chainid 1 and element != "H"') #this will calculate all residues against all residues. If you only want the interaciton between protein (chainid 0) and peptide(chainid 1), change the selection and the indices in the df.
protein_peptide2 = topology.select('chainid 0 or chainid 1 and element != "H"')

# 
print('Calculating Frequency of contacts')
contacts = ContactFrequency(traj, query=protein_peptide1, haystack=protein_peptide2, cutoff=0.45) # cutoff = size of the sphere in nm used to calculate the contacts in.
x = contacts.residue_contacts.sparse_matrix.toarray()
df = pd.DataFrame(x)
 
#NSD2 has 231 residues + 15 residues peptide = 246
df = df.iloc[: , :-4] # drop last 4 columns
df = df[:-4] # drop last 4 rows
if Peptide == 'H3K36':
    df.index = ['Y991', 'K992', 'H993', 'I994', 'K995', 'V996', 'N997', 'K998', 'P999', 'Y1000', 'G1001', 'K1002', 'V1003', 'Q1004', 'I1005', 'Y1006', 'T1007', 'A1008', 'D1009', 'I1010', 'S1011', 'E1012', 'I1013', 'P1014', 'K1015', 'C1016', 'N1017', 'C1018', 'K1019', 'P1020', 'T1021', 'D1022', 'E1023', 'N1024', 'P1025', 'C1026', 'G1027', 'F1028', 'D1029', 'S1030', 'E1031', 'C1032', 'L1033', 'N1034', 'R1035', 'M1036', 'L1037', 'M1038', 'F1039', 'E1040', 'H1041', 'P1043', 'Q1044', 'V1045', 'C1046', 'P1047', 'A1048', 'G1049', 'E1050', 'F1051', 'C1052', 'Q1053', 'N1054', 'Q1055', 'C1056', 'I1057', 'F1058', 'T1059', 'K1060', 'R1061', 'Q1062', 'Y1063', 'P1064', 'E1065', 'T1066', 'K1067', 'I1068', 'I1069', 'K1070', 'T1071', 'D1072', 'G1073', 'K1074', 'G1075', 'W1076', 'G1077', 'L1078', 'V1079', 'A1080', 'K1081', 'R1082', 'D1083', 'I1084', 'R1085', 'K1086', 'G1087', 'E1088', 'F1089', 'V1090', 'N1091', 'E1092', 'Y1093', 'V1094', 'G1095', 'E1096', 'L1097', 'I1098', 'D1099', 'E1100', 'E1101', 'E1102', 'C1103', 'M1104', 'A1105', 'R1106', 'I1107', 'K1108', 'H1109', 'A1110', 'H1111', 'E1112', 'N1113', 'D1114', 'I1115', 'T1116', 'H1117', 'F1118', 'Y1119', 'M1120', 'L1121', 'T1122', 'I1123', 'D1124', 'K1125', 'D1126', 'R1127', 'I1128', 'I1129', 'D1130', 'A1131', 'G1132', 'P1133', 'K1134', 'G1135', 'N1136', 'Y1137', 'S1138', 'R1139', 'F1140', 'M1141', 'N1142', 'H1143', 'S1144', 'C1145', 'Q1146', 'P1147', 'N1148', 'C1149', 'E1150', 'T1151', 'L1152', 'K1153', 'W1154', 'T1155', 'V1156', 'N1157', 'G1158', 'D1159', 'T1160', 'R1161', 'V1162', 'G1163', 'L1164', 'F1165', 'A1166', 'V1167', 'C1168', 'D1169', 'I1170', 'P1171', 'A1172', 'G1173', 'T1174', 'E1175', 'L1176', 'T1177', 'F1178', 'N1179', 'Y1180', 'N1181', 'L1182', 'D1183', 'C1184', 'L1185', 'G1186', 'N1187', 'E1188', 'K1189', 'T1190', 'V1191', 'C1192', 'R1193', 'C1194', 'G1195', 'A1196', 'S1197', 'N1198', 'C1199', 'S1200', 'G1201', 'F1202', 'L1203', 'G1204', 'D1205', 'R1206', 'P1207', 'K1208', 'T1209', 'S1210', 'T1211', 'T1212', 'L1213', 'S1214', 'S1215', 'E1216', 'E1217', 'K1218', 'G1219', 'K1220', 'K1221', 'A29', 'P30', 'A31', 'T32', 'G32', 'G33', 'V35', 'K36', 'K37', 'P38', 'H39', 'R40', 'Y41', 'R42', 'P43']
    df.columns = ['Y991', 'K992', 'H993', 'I994', 'K995', 'V996', 'N997', 'K998', 'P999', 'Y1000', 'G1001', 'K1002', 'V1003', 'Q1004', 'I1005', 'Y1006', 'T1007', 'A1008', 'D1009', 'I1010', 'S1011', 'E1012', 'I1013', 'P1014', 'K1015', 'C1016', 'N1017', 'C1018', 'K1019', 'P1020', 'T1021', 'D1022', 'E1023', 'N1024', 'P1025', 'C1026', 'G1027', 'F1028', 'D1029', 'S1030', 'E1031', 'C1032', 'L1033', 'N1034', 'R1035', 'M1036', 'L1037', 'M1038', 'F1039', 'E1040', 'H1041', 'P1043', 'Q1044', 'V1045', 'C1046', 'P1047', 'A1048', 'G1049', 'E1050', 'F1051', 'C1052', 'Q1053', 'N1054', 'Q1055', 'C1056', 'I1057', 'F1058', 'T1059', 'K1060', 'R1061', 'Q1062', 'Y1063', 'P1064', 'E1065', 'T1066', 'K1067', 'I1068', 'I1069', 'K1070', 'T1071', 'D1072', 'G1073', 'K1074', 'G1075', 'W1076', 'G1077', 'L1078', 'V1079', 'A1080', 'K1081', 'R1082', 'D1083', 'I1084', 'R1085', 'K1086', 'G1087', 'E1088', 'F1089', 'V1090', 'N1091', 'E1092', 'Y1093', 'V1094', 'G1095', 'E1096', 'L1097', 'I1098', 'D1099', 'E1100', 'E1101', 'E1102', 'C1103', 'M1104', 'A1105', 'R1106', 'I1107', 'K1108', 'H1109', 'A1110', 'H1111', 'E1112', 'N1113', 'D1114', 'I1115', 'T1116', 'H1117', 'F1118', 'Y1119', 'M1120', 'L1121', 'T1122', 'I1123', 'D1124', 'K1125', 'D1126', 'R1127', 'I1128', 'I1129', 'D1130', 'A1131', 'G1132', 'P1133', 'K1134', 'G1135', 'N1136', 'Y1137', 'S1138', 'R1139', 'F1140', 'M1141', 'N1142', 'H1143', 'S1144', 'C1145', 'Q1146', 'P1147', 'N1148', 'C1149', 'E1150', 'T1151', 'L1152', 'K1153', 'W1154', 'T1155', 'V1156', 'N1157', 'G1158', 'D1159', 'T1160', 'R1161', 'V1162', 'G1163', 'L1164', 'F1165', 'A1166', 'V1167', 'C1168', 'D1169', 'I1170', 'P1171', 'A1172', 'G1173', 'T1174', 'E1175', 'L1176', 'T1177', 'F1178', 'N1179', 'Y1180', 'N1181', 'L1182', 'D1183', 'C1184', 'L1185', 'G1186', 'N1187', 'E1188', 'K1189', 'T1190', 'V1191', 'C1192', 'R1193', 'C1194', 'G1195', 'A1196', 'S1197', 'N1198', 'C1199', 'S1200', 'G1201', 'F1202', 'L1203', 'G1204', 'D1205', 'R1206', 'P1207', 'K1208', 'T1209', 'S1210', 'T1211', 'T1212', 'L1213', 'S1214', 'S1215', 'E1216', 'E1217', 'K1218', 'G1219', 'K1220', 'K1221', 'A29', 'P30', 'A31', 'T32', 'G32', 'G33', 'V35', 'K36', 'K37', 'P38', 'H39', 'R40', 'Y41', 'R42', 'P43']
if Peptide == 'ssK36':
    df.index = ['Y991', 'K992', 'H993', 'I994', 'K995', 'V996', 'N997', 'K998', 'P999', 'Y1000', 'G1001', 'K1002', 'V1003', 'Q1004', 'I1005', 'Y1006', 'T1007', 'A1008', 'D1009', 'I1010', 'S1011', 'E1012', 'I1013', 'P1014', 'K1015', 'C1016', 'N1017', 'C1018', 'K1019', 'P1020', 'T1021', 'D1022', 'E1023', 'N1024', 'P1025', 'C1026', 'G1027', 'F1028', 'D1029', 'S1030', 'E1031', 'C1032', 'L1033', 'N1034', 'R1035', 'M1036', 'L1037', 'M1038', 'F1039', 'E1040', 'H1041', 'P1043', 'Q1044', 'V1045', 'C1046', 'P1047', 'A1048', 'G1049', 'E1050', 'F1051', 'C1052', 'Q1053', 'N1054', 'Q1055', 'C1056', 'I1057', 'F1058', 'T1059', 'K1060', 'R1061', 'Q1062', 'Y1063', 'P1064', 'E1065', 'T1066', 'K1067', 'I1068', 'I1069', 'K1070', 'T1071', 'D1072', 'G1073', 'K1074', 'G1075', 'W1076', 'G1077', 'L1078', 'V1079', 'A1080', 'K1081', 'R1082', 'D1083', 'I1084', 'R1085', 'K1086', 'G1087', 'E1088', 'F1089', 'V1090', 'N1091', 'E1092', 'Y1093', 'V1094', 'G1095', 'E1096', 'L1097', 'I1098', 'D1099', 'E1100', 'E1101', 'E1102', 'C1103', 'M1104', 'A1105', 'R1106', 'I1107', 'K1108', 'H1109', 'A1110', 'H1111', 'E1112', 'N1113', 'D1114', 'I1115', 'T1116', 'H1117', 'F1118', 'Y1119', 'M1120', 'L1121', 'T1122', 'I1123', 'D1124', 'K1125', 'D1126', 'R1127', 'I1128', 'I1129', 'D1130', 'A1131', 'G1132', 'P1133', 'K1134', 'G1135', 'N1136', 'Y1137', 'S1138', 'R1139', 'F1140', 'M1141', 'N1142', 'H1143', 'S1144', 'C1145', 'Q1146', 'P1147', 'N1148', 'C1149', 'E1150', 'T1151', 'L1152', 'K1153', 'W1154', 'T1155', 'V1156', 'N1157', 'G1158', 'D1159', 'T1160', 'R1161', 'V1162', 'G1163', 'L1164', 'F1165', 'A1166', 'V1167', 'C1168', 'D1169', 'I1170', 'P1171', 'A1172', 'G1173', 'T1174', 'E1175', 'L1176', 'T1177', 'F1178', 'N1179', 'Y1180', 'N1181', 'L1182', 'D1183', 'C1184', 'L1185', 'G1186', 'N1187', 'E1188', 'K1189', 'T1190', 'V1191', 'C1192', 'R1193', 'C1194', 'G1195', 'A1196', 'S1197', 'N1198', 'C1199', 'S1200', 'G1201', 'F1202', 'L1203', 'G1204', 'D1205', 'R1206', 'P1207', 'K1208', 'T1209', 'S1210', 'T1211', 'T1212', 'L1213', 'S1214', 'S1215', 'E1216', 'E1217', 'K1218', 'G1219', 'K1220', 'K1221', 'A29', 'P30', 'K31', 'T32', 'G32', 'G33', 'V35', 'K36', 'R37', 'P38', 'N39', 'N40', 'Y41', 'R42', 'P43']
    df.columns = ['Y991', 'K992', 'H993', 'I994', 'K995', 'V996', 'N997', 'K998', 'P999', 'Y1000', 'G1001', 'K1002', 'V1003', 'Q1004', 'I1005', 'Y1006', 'T1007', 'A1008', 'D1009', 'I1010', 'S1011', 'E1012', 'I1013', 'P1014', 'K1015', 'C1016', 'N1017', 'C1018', 'K1019', 'P1020', 'T1021', 'D1022', 'E1023', 'N1024', 'P1025', 'C1026', 'G1027', 'F1028', 'D1029', 'S1030', 'E1031', 'C1032', 'L1033', 'N1034', 'R1035', 'M1036', 'L1037', 'M1038', 'F1039', 'E1040', 'H1041', 'P1043', 'Q1044', 'V1045', 'C1046', 'P1047', 'A1048', 'G1049', 'E1050', 'F1051', 'C1052', 'Q1053', 'N1054', 'Q1055', 'C1056', 'I1057', 'F1058', 'T1059', 'K1060', 'R1061', 'Q1062', 'Y1063', 'P1064', 'E1065', 'T1066', 'K1067', 'I1068', 'I1069', 'K1070', 'T1071', 'D1072', 'G1073', 'K1074', 'G1075', 'W1076', 'G1077', 'L1078', 'V1079', 'A1080', 'K1081', 'R1082', 'D1083', 'I1084', 'R1085', 'K1086', 'G1087', 'E1088', 'F1089', 'V1090', 'N1091', 'E1092', 'Y1093', 'V1094', 'G1095', 'E1096', 'L1097', 'I1098', 'D1099', 'E1100', 'E1101', 'E1102', 'C1103', 'M1104', 'A1105', 'R1106', 'I1107', 'K1108', 'H1109', 'A1110', 'H1111', 'E1112', 'N1113', 'D1114', 'I1115', 'T1116', 'H1117', 'F1118', 'Y1119', 'M1120', 'L1121', 'T1122', 'I1123', 'D1124', 'K1125', 'D1126', 'R1127', 'I1128', 'I1129', 'D1130', 'A1131', 'G1132', 'P1133', 'K1134', 'G1135', 'N1136', 'Y1137', 'S1138', 'R1139', 'F1140', 'M1141', 'N1142', 'H1143', 'S1144', 'C1145', 'Q1146', 'P1147', 'N1148', 'C1149', 'E1150', 'T1151', 'L1152', 'K1153', 'W1154', 'T1155', 'V1156', 'N1157', 'G1158', 'D1159', 'T1160', 'R1161', 'V1162', 'G1163', 'L1164', 'F1165', 'A1166', 'V1167', 'C1168', 'D1169', 'I1170', 'P1171', 'A1172', 'G1173', 'T1174', 'E1175', 'L1176', 'T1177', 'F1178', 'N1179', 'Y1180', 'N1181', 'L1182', 'D1183', 'C1184', 'L1185', 'G1186', 'N1187', 'E1188', 'K1189', 'T1190', 'V1191', 'C1192', 'R1193', 'C1194', 'G1195', 'A1196', 'S1197', 'N1198', 'C1199', 'S1200', 'G1201', 'F1202', 'L1203', 'G1204', 'D1205', 'R1206', 'P1207', 'K1208', 'T1209', 'S1210', 'T1211', 'T1212', 'L1213', 'S1214', 'S1215', 'E1216', 'E1217', 'K1218', 'G1219', 'K1220', 'K1221', 'A29', 'P30', 'K31', 'T32', 'G32', 'G33', 'V35', 'K36', 'R37', 'P38', 'N39', 'N40', 'Y41', 'R42', 'P43']

# create new folder and safe pickle, excel sheet
new_folder = 'contacts_{}_{}'.format(Protein, Peptide)
os.system('mkdir {0}'.format(new_folder))
df.to_pickle(open(new_folder + '/' +  'df_contacts_{}_{}.pkl'.format(Protein, Peptide), 'wb'))
df.to_excel(open(new_folder + '/' +  'contacts_{}_{}.xlsx'.format(Protein, Peptide), 'wb'))
