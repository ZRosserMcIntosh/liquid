// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
//
// Minimal JSON (RFC 8259) DOM parser for liquid, surface, vessel and projectile definitions.
// Dependency-free so the core builds unchanged for WebAssembly. Errors carry
// line and column. Duplicate object keys are rejected.
#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace lqd::json {

class Value {
public:
    enum class Type { Null, Bool, Number, String, Array, Object };

    Value() = default;
    static Value number(double d) { Value v; v.type_ = Type::Number; v.num_ = d; return v; }
    static Value boolean(bool b) { Value v; v.type_ = Type::Bool; v.bool_ = b; return v; }
    static Value string(std::string s) { Value v; v.type_ = Type::String; v.str_ = std::move(s); return v; }
    static Value array() { Value v; v.type_ = Type::Array; return v; }
    static Value object() { Value v; v.type_ = Type::Object; return v; }

    Type type() const { return type_; }
    bool isNull() const { return type_ == Type::Null; }
    bool isNumber() const { return type_ == Type::Number; }
    bool isBool() const { return type_ == Type::Bool; }
    bool isString() const { return type_ == Type::String; }
    bool isArray() const { return type_ == Type::Array; }
    bool isObject() const { return type_ == Type::Object; }

    double asNumber() const { return num_; }
    bool asBool() const { return bool_; }
    const std::string& asString() const { return str_; }
    const std::vector<Value>& items() const { return arr_; }
    std::vector<Value>& items() { return arr_; }
    // Object members in source order.
    const std::vector<std::pair<std::string, Value>>& members() const { return obj_; }
    std::vector<std::pair<std::string, Value>>& members() { return obj_; }
    const Value* find(const std::string& key) const;

private:
    Type type_ = Type::Null;
    double num_ = 0;
    bool bool_ = false;
    std::string str_;
    std::vector<Value> arr_;
    std::vector<std::pair<std::string, Value>> obj_;
};

struct ParseResult {
    Value value;
    std::string error;  // empty on success; "line L, column C: message"
    bool ok() const { return error.empty(); }
};

ParseResult parse(const std::string& text);

}  // namespace lqd::json
