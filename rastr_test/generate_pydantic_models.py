import json
from datamodel_code_generator import InputFileType, generate
# 1. Импортируем перечисление Formatter и PythonVersion
from datamodel_code_generator.format import Formatter
from datamodel_code_generator.types import PythonVersion
import re
import astra


mode_from_rastr = True


def transliterate(text: str) -> str:
    """Переводит кириллицу в латиницу и делает строку валидным ID для Python."""
    cyrillic = "абвгдеёжзийклмнопрстуфхцчшщъыьэюя"
    latin = "abvgdeezzijklmnoprstufhcshsh''eua"
    tr_map = str.maketrans(cyrillic, latin)

    text = text.lower().translate(tr_map)
    # Заменяем пробелы, дефисы и спецсимволы на нижнее подчеркивание
    text = re.sub(r"[^a-z0-9_]", "_", text)
    # Если строка начинается с цифры, добавляем префикс val_
    if text and text[0].isdigit():
        text = f"val_{text}"
    return text or "empty"

def generate_model() -> json:
    rastr = astra.Rastr()
    json_schema = {
    "$schema": "http://json-schema.org",
    "title": "RastrModel",
    "type": "object",
    "properties": {
        "version": {
                "type": "string",
                "default": "1.0.0",
                "const": "1.0.0",
                "description": "Версия модели"
            }
        },
    "$defs": {},
    "additionalProperties": False,
    "required": ["version"]
 }
    type_map = {
        astra.PropType.INT: "integer", 
        astra.PropType.DOUBLE: "number",
        astra.PropType.STRING: "string",
        astra.PropType.BOOL: "boolean",
        astra.PropType.ENUM: "string"
        }

    rastr.new_file("d:/documents/rastrwin3\shablon\режим.rg2")
    for table in rastr.tables():
        table_keys = table.key.split(",")
        table_data = {
            "type": "object",
            "additionalProperties": False
        }
        fields_data = {}
        required = []
        for field in table.columns():

            if field.name.startswith("_"):
                continue

            add_to_required = field.name in table_keys

            read_only_field = len(field.property(astra.FieldProperties.EXPRESSION)) > 0

            if mode_from_rastr:
                # отдаем вычисляемые поля
                add_to_required = True
            else:
                # не принимаем вычисляемые поля
                if read_only_field:
                    continue


            if field.type in type_map:
                field_data = {
                    "description" : field.property(astra.FieldProperties.DESCRIPTION),
                    "title" : field.property(astra.FieldProperties.TITLE),
                }
                if field.type == astra.PropType.ENUM:
                    enum_items = field.property(astra.FieldProperties.NAMEREF).split("|")
                    enum_name = f"{table.name}{field.name}Enum"
                
                    json_schema["$defs"][enum_name] = {
                        "type": "string",
                        "title": enum_name,
                        "enum": enum_items,
                        "x-enum-varnames": [transliterate(x) for x in enum_items]
                    }
                    field_data["$ref"] = f"#/$defs/{enum_name}"
                else:
                    field_data["type"] = type_map[field.type]


                fields_data[field.name] = field_data
                add_to_required |= field.type == astra.PropType.BOOL
                add_to_required |= len(field.property(astra.FieldProperties.NAMEREF)) and field.type == astra.PropType.INT
                if add_to_required:
                    required.append(field.name)
                        
        table_data["properties"] = fields_data
        table_data["required"] = required
        json_schema["properties"][table.name] = {"type": "array", "items": table_data}



    python_code = generate(
        json.dumps(json_schema),
        input_file_type=InputFileType.JsonSchema,
        output_model_type="pydantic_v2.BaseModel",
        target_python_version=PythonVersion.PY_310,
        formatters=[Formatter.BLACK, Formatter.ISORT]
    )
    with open("models.py", "w", encoding="utf-8") as f:

        f.write(python_code)

    return {}

generate_model()
    
