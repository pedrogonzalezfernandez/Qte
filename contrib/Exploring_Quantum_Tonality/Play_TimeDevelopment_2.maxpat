{
	"patcher" : 	{
		"fileversion" : 1,
		"appversion" : 		{
			"major" : 9,
			"minor" : 0,
			"revision" : 8,
			"architecture" : "x64",
			"modernui" : 1
		}
,
		"classnamespace" : "box",
		"rect" : [ 166.0, 154.0, 953.0, 405.0 ],
		"gridsize" : [ 15.0, 15.0 ],
		"boxes" : [ 			{
				"box" : 				{
					"id" : "obj-55",
					"maxclass" : "newobj",
					"numinlets" : 1,
					"numoutlets" : 1,
					"outlettype" : [ "signal" ],
					"patching_rect" : [ 228.0, 302.0, 110.0, 22.0 ],
					"text" : "receive~ master_R"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-54",
					"maxclass" : "newobj",
					"numinlets" : 1,
					"numoutlets" : 1,
					"outlettype" : [ "signal" ],
					"patching_rect" : [ 43.0, 302.0, 108.0, 22.0 ],
					"text" : "receive~ master_L"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-53",
					"maxclass" : "ezdac~",
					"numinlets" : 2,
					"numoutlets" : 0,
					"patching_rect" : [ 160.0, 342.0, 45.0, 45.0 ]
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-44",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 0,
					"patching_rect" : [ 507.0, 234.0, 78.0, 22.0 ],
					"text" : "PlayVoice 16"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-45",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 0,
					"patching_rect" : [ 426.0, 234.0, 78.0, 22.0 ],
					"text" : "PlayVoice 15"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-46",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 0,
					"patching_rect" : [ 348.0, 234.0, 78.0, 22.0 ],
					"text" : "PlayVoice 14"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-47",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 0,
					"patching_rect" : [ 268.0, 234.0, 78.0, 22.0 ],
					"text" : "PlayVoice 13"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-48",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 0,
					"patching_rect" : [ 188.0, 234.0, 78.0, 22.0 ],
					"text" : "PlayVoice 12"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-49",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 0,
					"patching_rect" : [ 103.0, 234.0, 77.0, 22.0 ],
					"text" : "PlayVoice 11"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-50",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 0,
					"patching_rect" : [ 17.0, 234.0, 78.0, 22.0 ],
					"text" : "PlayVoice 10"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-51",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 0,
					"patching_rect" : [ 589.0, 234.0, 78.0, 22.0 ],
					"text" : "PlayVoice 17"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-43",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 0,
					"patching_rect" : [ 564.0, 89.0, 71.0, 22.0 ],
					"text" : "PlayVoice 8"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-42",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 0,
					"patching_rect" : [ 488.0, 89.0, 71.0, 22.0 ],
					"text" : "PlayVoice 7"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-40",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 0,
					"patching_rect" : [ 412.0, 89.0, 71.0, 22.0 ],
					"text" : "PlayVoice 6"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-39",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 0,
					"patching_rect" : [ 331.0, 89.0, 71.0, 22.0 ],
					"text" : "PlayVoice 5"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-38",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 0,
					"patching_rect" : [ 251.0, 89.0, 71.0, 22.0 ],
					"text" : "PlayVoice 4"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-5",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 0,
					"patching_rect" : [ 172.0, 89.0, 71.0, 22.0 ],
					"text" : "PlayVoice 3"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-4",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 0,
					"patching_rect" : [ 94.0, 89.0, 71.0, 22.0 ],
					"text" : "PlayVoice 2"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-3",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 0,
					"patching_rect" : [ 17.0, 89.0, 71.0, 22.0 ],
					"text" : "PlayVoice 1"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-41",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 0,
					"patching_rect" : [ 642.0, 89.0, 71.0, 22.0 ],
					"text" : "PlayVoice 9"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-24",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 533.0, 180.0, 69.0, 22.0 ],
					"text" : "r 17_phase"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-25",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 518.0, 156.0, 59.0, 22.0 ],
					"text" : "r 17_mag"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-26",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 461.0, 180.0, 69.0, 22.0 ],
					"text" : "r 16_phase"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-27",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 446.0, 156.0, 59.0, 22.0 ],
					"text" : "r 16_mag"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-28",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 392.0, 180.0, 69.0, 22.0 ],
					"text" : "r 15_phase"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-29",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 377.0, 156.0, 59.0, 22.0 ],
					"text" : "r 15_mag"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-30",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 325.0, 180.0, 69.0, 22.0 ],
					"text" : "r 14_phase"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-31",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 310.0, 156.0, 59.0, 22.0 ],
					"text" : "r 14_mag"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-32",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 253.0, 180.0, 69.0, 22.0 ],
					"text" : "r 13_phase"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-33",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 238.0, 156.0, 59.0, 22.0 ],
					"text" : "r 13_mag"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-34",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 182.0, 180.0, 69.0, 22.0 ],
					"text" : "r 12_phase"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-35",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 167.0, 156.0, 59.0, 22.0 ],
					"text" : "r 12_mag"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-36",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 110.0, 180.0, 68.0, 22.0 ],
					"text" : "r 11_phase"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-37",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 95.0, 156.0, 58.0, 22.0 ],
					"text" : "r 11_mag"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-14",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 33.0, 180.0, 69.0, 22.0 ],
					"text" : "r 10_phase"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-15",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 18.0, 156.0, 59.0, 22.0 ],
					"text" : "r 10_mag"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-16",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 594.0, 46.0, 62.0, 22.0 ],
					"text" : "r 9_phase"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-17",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 579.0, 22.0, 53.0, 22.0 ],
					"text" : "r 9_mag"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-18",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 522.0, 46.0, 62.0, 22.0 ],
					"text" : "r 8_phase"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-19",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 507.0, 22.0, 53.0, 22.0 ],
					"text" : "r 8_mag"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-20",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 451.0, 46.0, 62.0, 22.0 ],
					"text" : "r 7_phase"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-21",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 436.0, 22.0, 53.0, 22.0 ],
					"text" : "r 7_mag"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-22",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 379.0, 46.0, 62.0, 22.0 ],
					"text" : "r 6_phase"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-23",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 364.0, 22.0, 53.0, 22.0 ],
					"text" : "r 6_mag"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-12",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 310.0, 46.0, 62.0, 22.0 ],
					"text" : "r 5_phase"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-13",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 295.0, 22.0, 53.0, 22.0 ],
					"text" : "r 5_mag"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-10",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 243.0, 46.0, 62.0, 22.0 ],
					"text" : "r 4_phase"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-11",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 228.0, 22.0, 53.0, 22.0 ],
					"text" : "r 4_mag"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-8",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 171.0, 46.0, 62.0, 22.0 ],
					"text" : "r 3_phase"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-9",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 156.0, 22.0, 53.0, 22.0 ],
					"text" : "r 3_mag"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-6",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 100.0, 46.0, 62.0, 22.0 ],
					"text" : "r 2_phase"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-7",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 85.0, 22.0, 53.0, 22.0 ],
					"text" : "r 2_mag"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-1",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 32.0, 46.0, 62.0, 22.0 ],
					"text" : "r 1_phase"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-2",
					"maxclass" : "newobj",
					"numinlets" : 0,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 24.0, 22.0, 53.0, 22.0 ],
					"text" : "r 1_mag"
				}

			}
 ],
		"lines" : [ 			{
				"patchline" : 				{
					"destination" : [ "obj-3", 1 ],
					"source" : [ "obj-1", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-38", 1 ],
					"source" : [ "obj-10", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-38", 0 ],
					"source" : [ "obj-11", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-39", 1 ],
					"source" : [ "obj-12", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-39", 0 ],
					"source" : [ "obj-13", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-50", 1 ],
					"source" : [ "obj-14", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-50", 0 ],
					"source" : [ "obj-15", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-41", 1 ],
					"source" : [ "obj-16", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-41", 0 ],
					"source" : [ "obj-17", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-43", 1 ],
					"source" : [ "obj-18", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-43", 0 ],
					"source" : [ "obj-19", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-3", 0 ],
					"source" : [ "obj-2", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-42", 1 ],
					"source" : [ "obj-20", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-42", 0 ],
					"source" : [ "obj-21", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-40", 1 ],
					"source" : [ "obj-22", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-40", 0 ],
					"source" : [ "obj-23", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-51", 1 ],
					"source" : [ "obj-24", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-51", 0 ],
					"source" : [ "obj-25", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-44", 1 ],
					"source" : [ "obj-26", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-44", 0 ],
					"source" : [ "obj-27", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-45", 1 ],
					"source" : [ "obj-28", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-45", 0 ],
					"source" : [ "obj-29", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-46", 1 ],
					"source" : [ "obj-30", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-46", 0 ],
					"source" : [ "obj-31", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-47", 1 ],
					"source" : [ "obj-32", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-47", 0 ],
					"source" : [ "obj-33", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-48", 1 ],
					"source" : [ "obj-34", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-48", 0 ],
					"source" : [ "obj-35", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-49", 1 ],
					"source" : [ "obj-36", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-49", 0 ],
					"source" : [ "obj-37", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-53", 0 ],
					"source" : [ "obj-54", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-53", 1 ],
					"source" : [ "obj-55", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-4", 1 ],
					"source" : [ "obj-6", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-4", 0 ],
					"source" : [ "obj-7", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-5", 1 ],
					"source" : [ "obj-8", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-5", 0 ],
					"source" : [ "obj-9", 0 ]
				}

			}
 ],
		"parameters" : 		{
			"obj-38::obj-2" : [ "number[12]", "number[1]", 0 ],
			"obj-38::obj-35" : [ "number[16]", "number[4]", 0 ],
			"obj-38::obj-36" : [ "number[15]", "number[3]", 0 ],
			"obj-38::obj-4" : [ "number[11]", "number[2]", 0 ],
			"obj-38::obj-60" : [ "number[13]", "number[6]", 0 ],
			"obj-38::obj-61" : [ "number[14]", "number[5]", 0 ],
			"obj-38::obj-64" : [ "number[10]", "number[8]", 0 ],
			"obj-38::obj-78" : [ "number[189]", "number[9]", 0 ],
			"obj-39::obj-2" : [ "number[183]", "number[1]", 0 ],
			"obj-39::obj-35" : [ "number[185]", "number[4]", 0 ],
			"obj-39::obj-36" : [ "number[182]", "number[3]", 0 ],
			"obj-39::obj-4" : [ "number[186]", "number[2]", 0 ],
			"obj-39::obj-60" : [ "number[184]", "number[6]", 0 ],
			"obj-39::obj-61" : [ "number[181]", "number[5]", 0 ],
			"obj-39::obj-64" : [ "number[188]", "number[8]", 0 ],
			"obj-39::obj-78" : [ "number[187]", "number[9]", 0 ],
			"obj-3::obj-2" : [ "number[1]", "number[1]", 0 ],
			"obj-3::obj-35" : [ "number[4]", "number[4]", 0 ],
			"obj-3::obj-36" : [ "number[3]", "number[3]", 0 ],
			"obj-3::obj-4" : [ "number[2]", "number[2]", 0 ],
			"obj-3::obj-60" : [ "number[6]", "number[6]", 0 ],
			"obj-3::obj-61" : [ "number[5]", "number[5]", 0 ],
			"obj-3::obj-64" : [ "number[8]", "number[8]", 0 ],
			"obj-3::obj-78" : [ "number[9]", "number[9]", 0 ],
			"obj-40::obj-2" : [ "number[180]", "number[1]", 0 ],
			"obj-40::obj-35" : [ "number[174]", "number[4]", 0 ],
			"obj-40::obj-36" : [ "number[179]", "number[3]", 0 ],
			"obj-40::obj-4" : [ "number[175]", "number[2]", 0 ],
			"obj-40::obj-60" : [ "number[173]", "number[6]", 0 ],
			"obj-40::obj-61" : [ "number[178]", "number[5]", 0 ],
			"obj-40::obj-64" : [ "number[177]", "number[8]", 0 ],
			"obj-40::obj-78" : [ "number[176]", "number[9]", 0 ],
			"obj-41::obj-2" : [ "number[34]", "number[1]", 0 ],
			"obj-41::obj-35" : [ "number[40]", "number[4]", 0 ],
			"obj-41::obj-36" : [ "number[38]", "number[3]", 0 ],
			"obj-41::obj-4" : [ "number[33]", "number[2]", 0 ],
			"obj-41::obj-60" : [ "number[39]", "number[6]", 0 ],
			"obj-41::obj-61" : [ "number[37]", "number[5]", 0 ],
			"obj-41::obj-64" : [ "number[36]", "number[8]", 0 ],
			"obj-41::obj-78" : [ "number[35]", "number[9]", 0 ],
			"obj-42::obj-2" : [ "number[169]", "number[1]", 0 ],
			"obj-42::obj-35" : [ "number[171]", "number[4]", 0 ],
			"obj-42::obj-36" : [ "number[170]", "number[3]", 0 ],
			"obj-42::obj-4" : [ "number[168]", "number[2]", 0 ],
			"obj-42::obj-60" : [ "number[167]", "number[6]", 0 ],
			"obj-42::obj-61" : [ "number[172]", "number[5]", 0 ],
			"obj-42::obj-64" : [ "number[166]", "number[8]", 0 ],
			"obj-42::obj-78" : [ "number[165]", "number[9]", 0 ],
			"obj-43::obj-2" : [ "number[158]", "number[1]", 0 ],
			"obj-43::obj-35" : [ "number[161]", "number[4]", 0 ],
			"obj-43::obj-36" : [ "number[159]", "number[3]", 0 ],
			"obj-43::obj-4" : [ "number[157]", "number[2]", 0 ],
			"obj-43::obj-60" : [ "number[160]", "number[6]", 0 ],
			"obj-43::obj-61" : [ "number[162]", "number[5]", 0 ],
			"obj-43::obj-64" : [ "number[163]", "number[8]", 0 ],
			"obj-43::obj-78" : [ "number[164]", "number[9]", 0 ],
			"obj-44::obj-2" : [ "number[94]", "number[1]", 0 ],
			"obj-44::obj-35" : [ "number[93]", "number[4]", 0 ],
			"obj-44::obj-36" : [ "number[98]", "number[3]", 0 ],
			"obj-44::obj-4" : [ "number[96]", "number[2]", 0 ],
			"obj-44::obj-60" : [ "number[97]", "number[6]", 0 ],
			"obj-44::obj-61" : [ "number[95]", "number[5]", 0 ],
			"obj-44::obj-64" : [ "number[99]", "number[8]", 0 ],
			"obj-44::obj-78" : [ "number[100]", "number[9]", 0 ],
			"obj-45::obj-2" : [ "number[104]", "number[1]", 0 ],
			"obj-45::obj-35" : [ "number[101]", "number[4]", 0 ],
			"obj-45::obj-36" : [ "number[103]", "number[3]", 0 ],
			"obj-45::obj-4" : [ "number[106]", "number[2]", 0 ],
			"obj-45::obj-60" : [ "number[107]", "number[6]", 0 ],
			"obj-45::obj-61" : [ "number[108]", "number[5]", 0 ],
			"obj-45::obj-64" : [ "number[102]", "number[8]", 0 ],
			"obj-45::obj-78" : [ "number[105]", "number[9]", 0 ],
			"obj-46::obj-2" : [ "number[111]", "number[1]", 0 ],
			"obj-46::obj-35" : [ "number[115]", "number[4]", 0 ],
			"obj-46::obj-36" : [ "number[113]", "number[3]", 0 ],
			"obj-46::obj-4" : [ "number[116]", "number[2]", 0 ],
			"obj-46::obj-60" : [ "number[114]", "number[6]", 0 ],
			"obj-46::obj-61" : [ "number[112]", "number[5]", 0 ],
			"obj-46::obj-64" : [ "number[109]", "number[8]", 0 ],
			"obj-46::obj-78" : [ "number[110]", "number[9]", 0 ],
			"obj-47::obj-2" : [ "number[122]", "number[1]", 0 ],
			"obj-47::obj-35" : [ "number[121]", "number[4]", 0 ],
			"obj-47::obj-36" : [ "number[120]", "number[3]", 0 ],
			"obj-47::obj-4" : [ "number[123]", "number[2]", 0 ],
			"obj-47::obj-60" : [ "number[124]", "number[6]", 0 ],
			"obj-47::obj-61" : [ "number[119]", "number[5]", 0 ],
			"obj-47::obj-64" : [ "number[117]", "number[8]", 0 ],
			"obj-47::obj-78" : [ "number[118]", "number[9]", 0 ],
			"obj-48::obj-2" : [ "number[126]", "number[1]", 0 ],
			"obj-48::obj-35" : [ "number[130]", "number[4]", 0 ],
			"obj-48::obj-36" : [ "number[129]", "number[3]", 0 ],
			"obj-48::obj-4" : [ "number[128]", "number[2]", 0 ],
			"obj-48::obj-60" : [ "number[132]", "number[6]", 0 ],
			"obj-48::obj-61" : [ "number[131]", "number[5]", 0 ],
			"obj-48::obj-64" : [ "number[127]", "number[8]", 0 ],
			"obj-48::obj-78" : [ "number[125]", "number[9]", 0 ],
			"obj-49::obj-2" : [ "number[134]", "number[1]", 0 ],
			"obj-49::obj-35" : [ "number[135]", "number[4]", 0 ],
			"obj-49::obj-36" : [ "number[133]", "number[3]", 0 ],
			"obj-49::obj-4" : [ "number[136]", "number[2]", 0 ],
			"obj-49::obj-60" : [ "number[140]", "number[6]", 0 ],
			"obj-49::obj-61" : [ "number[137]", "number[5]", 0 ],
			"obj-49::obj-64" : [ "number[139]", "number[8]", 0 ],
			"obj-49::obj-78" : [ "number[138]", "number[9]", 0 ],
			"obj-4::obj-2" : [ "number[27]", "number[1]", 0 ],
			"obj-4::obj-35" : [ "number[32]", "number[4]", 0 ],
			"obj-4::obj-36" : [ "number[30]", "number[3]", 0 ],
			"obj-4::obj-4" : [ "number[28]", "number[2]", 0 ],
			"obj-4::obj-60" : [ "number[31]", "number[6]", 0 ],
			"obj-4::obj-61" : [ "number[25]", "number[5]", 0 ],
			"obj-4::obj-64" : [ "number[26]", "number[8]", 0 ],
			"obj-4::obj-78" : [ "number[29]", "number[9]", 0 ],
			"obj-50::obj-2" : [ "number[146]", "number[1]", 0 ],
			"obj-50::obj-35" : [ "number[145]", "number[4]", 0 ],
			"obj-50::obj-36" : [ "number[148]", "number[3]", 0 ],
			"obj-50::obj-4" : [ "number[147]", "number[2]", 0 ],
			"obj-50::obj-60" : [ "number[144]", "number[6]", 0 ],
			"obj-50::obj-61" : [ "number[143]", "number[5]", 0 ],
			"obj-50::obj-64" : [ "number[142]", "number[8]", 0 ],
			"obj-50::obj-78" : [ "number[141]", "number[9]", 0 ],
			"obj-51::obj-2" : [ "number[154]", "number[1]", 0 ],
			"obj-51::obj-35" : [ "number[155]", "number[4]", 0 ],
			"obj-51::obj-36" : [ "number[153]", "number[3]", 0 ],
			"obj-51::obj-4" : [ "number[156]", "number[2]", 0 ],
			"obj-51::obj-60" : [ "number[152]", "number[6]", 0 ],
			"obj-51::obj-61" : [ "number[149]", "number[5]", 0 ],
			"obj-51::obj-64" : [ "number[151]", "number[8]", 0 ],
			"obj-51::obj-78" : [ "number[150]", "number[9]", 0 ],
			"obj-5::obj-2" : [ "number[24]", "number[1]", 0 ],
			"obj-5::obj-35" : [ "number[23]", "number[4]", 0 ],
			"obj-5::obj-36" : [ "number[21]", "number[3]", 0 ],
			"obj-5::obj-4" : [ "number[18]", "number[2]", 0 ],
			"obj-5::obj-60" : [ "number[20]", "number[6]", 0 ],
			"obj-5::obj-61" : [ "number[17]", "number[5]", 0 ],
			"obj-5::obj-64" : [ "number[22]", "number[8]", 0 ],
			"obj-5::obj-78" : [ "number[19]", "number[9]", 0 ],
			"parameterbanks" : 			{
				"0" : 				{
					"index" : 0,
					"name" : "",
					"parameters" : [ "-", "-", "-", "-", "-", "-", "-", "-" ]
				}

			}
,
			"inherited_shortname" : 1
		}
,
		"dependency_cache" : [ 			{
				"name" : "PlayVoice.maxpat",
				"bootpath" : "~/Documents/Arbeit 2026/Palermo Proceedings/Final Submission/Programms",
				"patcherrelativepath" : ".",
				"type" : "JSON",
				"implicit" : 1
			}
 ],
		"autosave" : 0
	}

}
